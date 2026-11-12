//! `ul_realignment_gfa` (`gfa_ut.cpp` 17598) plus the `renew_g` /
//! `gradually_renew_g` drivers (`Overlaps.cpp` 39262 / 39297).

use crate::assembly::asg::Asg;
use crate::assembly::bubble::BubbleLabel;
use crate::assembly::hit::{CoverageCut, ReadToUnitig};
use crate::correction::ec_driver::PafList;
use crate::assembly::ul_store::AllUl;
use crate::command_line::options::AssemblyOptions;
use crate::reads::ReadStore;

/// The `ulg_opt_t` fields plus the round bookkeeping of the renew loop.
pub struct UlRenewOptions {
    /// `clean_round`.
    pub clean_round: i64,
    /// `min_ovlp_drop_ratio`.
    pub min_ovlp_drop_ratio: f64,
    /// `max_ovlp_drop_ratio`.
    pub max_ovlp_drop_ratio: f64,
    /// `asm_opt.max_short_tip`.
    pub max_tip: i64,
    /// `asm_opt.max_short_ul_tip`.
    pub max_ul_tip: i64,
    /// `ha_opt_triobin(&asm_opt)`.
    pub is_trio: bool,
}

/// `ul_realignment_gfa` (`gfa_ut.cpp` 17598): one renewal round — rebuild
/// the unitig graph against the UL alignments, regenerate the integer
/// graph, thread it and swap the read graph for the renewed one.
///
/// `deep_clean` (`deep_graph_clean`) and the trailing `hic_clean_adv` rely
/// on the hic.cpp bubble machinery, which is not ported; both stay no-ops
/// exactly like the non-`ar` stub in `pipeline.rs`.
#[allow(clippy::too_many_arguments)]
pub fn ul_realignment_gfa(
    options: &mut AssemblyOptions,
    reads: &mut ReadStore,
    sg: &mut Asg,
    sources: &mut Vec<PafList>,
    reverse_sources: &mut Vec<PafList>,
    coverage_cut: &mut Vec<CoverageCut>,
    read_to_unitig: &mut ReadToUnitig,
    bubble_labels: &mut BubbleLabel,
    ul_store: &mut AllUl,
    opts: &UlRenewOptions,
    bin_file: &str,
    _free_uld: bool,
    is_bridg: bool,
    deep_clean: bool,
) {
    for seq in sg.sequences.iter_mut() {
        if !seq.deleted {
            seq.coverage_class = crate::assembly::hit::PRIMARY_LABEL;
        }
    }
    crate::assembly::unitig::hic_clean(sg);
    let _ = deep_clean;

    // `init_ug = ul_realignment(uopt, sg, 0, bin_file)`
    let mut init_ug = crate::assembly::ul_recall::ul_realignment(
        options,
        reads,
        sources,
        coverage_cut,
        sg,
        bin_file,
        0,
        options.minimizer_sample_distance,
        options.minimizer_rewindow_size,
        options.maximum_hang_length,
        options.maximum_hang_rate,
        options.minimum_overlap_length,
    );

    // `filter_sg_by_ug(sg, init_ug, uopt)`
    crate::assembly::ul_renew_gfa::filter_sg_by_ug(
        sg,
        &mut init_ug,
        sources,
        options.maximum_hang_length,
        options.maximum_hang_rate,
        options.minimum_overlap_length,
    );

    // `bub = gen_bubble_chain(sg, init_ug, uopt, &r_het, polyploidy > 2)`
    let (mut bub, _r_het) = crate::assembly::ul_renew_gfa::gen_bubble_chain(
        sg,
        &mut init_ug,
        sources,
        reverse_sources,
        read_to_unitig,
        options,
        coverage_cut,
        Some(bubble_labels),
        options.polyploidy > 2,
    );

    // `init_ul_resolve_t`: the private `l1_ug` copy plus the string index.
    let mut l1_ug = crate::assembly::unitig::copy_untig_graph(&init_ug);
    let mut str_idx =
        crate::assembly::ul_integer::init_ul_str_idx(&ul_store.reads, l1_ug.unitigs.len());

    // `ul_re_correct(uidx, integer_correct_round)`: with the default zero
    // rounds only the `shrink_uls` tail runs.
    let hifi_read_count = reads.total_reads as usize;
    crate::assembly::ul_renew_ng::shrink_uls(
        ul_store,
        &mut str_idx,
        &mut l1_ug,
        reads,
        hifi_read_count,
    );

    // `init_ulg_opt_t(&uu, uopt, clean_round, ...)`
    let ulg_options = crate::assembly::ulg_clean::UlgOptions {
        clean_round: opts.clean_round,
        min_path_drop_ratio: f64::from(options.minimum_path_drop_rate),
        max_path_drop_ratio: f64::from(options.maximum_path_drop_rate),
        min_ovlp_drop_ratio: opts.min_ovlp_drop_ratio,
        max_ovlp_drop_ratio: opts.max_ovlp_drop_ratio,
        max_tip: opts.max_tip,
        max_tip_hifi: opts.max_ul_tip,
        is_trio: opts.is_trio,
    };

    // `gen_ul2ul(uidx, uopt, &uu, 0, is_bridg)`
    let hit_lists: Vec<Vec<crate::overlap::paf::PafHit>> =
        sources.iter().map(|list| list.hits.clone()).collect();
    let mut z = crate::assembly::ul2ul::gen_ul2ul(
        ul_store,
        &mut l1_ug,
        sg,
        &mut str_idx,
        &bub,
        &hit_lists,
        reads,
        i64::from(options.minimum_overlap_length),
        options.maximum_hang_length,
        options.maximum_hang_rate,
        options.gap_fuzz,
        bubble_labels,
        &ulg_options,
        false,
        is_bridg,
    );

    // `ng = renew_ng(uidx->uovl.h_usg, init_ug, sg, uopt, cov, ruIndex, 16, bin_file)`
    let h_usg = std::mem::take(&mut z.h_usg);
    let ng = crate::assembly::ul_renew_ng::renew_ng(
        &h_usg,
        &init_ug,
        sg,
        sources,
        reverse_sources,
        coverage_cut,
        read_to_unitig,
        reads,
        ul_store,
        options.maximum_hang_length,
        options.maximum_hang_rate,
        i64::from(options.minimum_overlap_length),
        16,
        &options.output_file_name,
        bin_file,
        bubble_labels,
    );

    // `(*(ropt->sg)) = ng`
    *sg = ng;

    // `renew_R_to_U(ng, src, r_src, n_read, cov, ruIndex, max_hang, min_ovlp)`
    crate::assembly::ul_renew_ng::renew_r_to_u(
        sources,
        &mut [],
        reverse_sources,
        &mut [],
        coverage_cut,
        read_to_unitig,
        options.maximum_hang_length,
        options.minimum_overlap_length,
        reads,
    );
    // `post_rescue` runs in the caller.
}

/// `gradually_renew_g` (`Overlaps.cpp` 39297): the renewal rounds with the
/// `ug_ext_gfa` sweeps in between.
#[allow(clippy::too_many_arguments)]
pub fn gradually_renew_g(
    options: &mut AssemblyOptions,
    reads: &mut ReadStore,
    sg: &mut Asg,
    sources: &mut Vec<PafList>,
    reverse_sources: &mut Vec<PafList>,
    coverage_cut: &mut Vec<CoverageCut>,
    read_to_unitig: &mut ReadToUnitig,
    bubble_labels: &mut BubbleLabel,
    ul_store: &mut AllUl,
    opts: &UlRenewOptions,
    ul_aln_round: i32,
    maximum_length: u32,
) {
    let bin_file = "re";
    ul_realignment_gfa(
        options,
        reads,
        sg,
        sources,
        reverse_sources,
        coverage_cut,
        read_to_unitig,
        bubble_labels,
        ul_store,
        opts,
        bin_file,
        ul_aln_round <= 1,
        true,
        options.polyploidy <= 2,
    );
    // `gen_ug_opt_t(..., (max_short_tip*2), 0.15, 3, 0.05, 0.9, ...)` then
    // `ug_ext_gfa(uopt, *sg, ug_ext_len)`
    crate::assembly::unitig::ug_ext_gfa(
        sg,
        sources,
        read_to_unitig,
        (opts.max_tip * 2) as u32,
        options.maximum_hang_length,
        options.minimum_overlap_length,
        options.maximum_hang_rate,
        maximum_length,
    );
    // `hic_clean_adv(*sg, uopt)` — not ported (see the module comment).

    for k in 1..ul_aln_round {
        let id = format!("{}{}", bin_file, k);
        ul_realignment_gfa(
            options,
            reads,
            sg,
            sources,
            reverse_sources,
            coverage_cut,
            read_to_unitig,
            bubble_labels,
            ul_store,
            opts,
            &id,
            (k + 1) == ul_aln_round,
            false,
            false,
        );
        crate::assembly::unitig::ug_ext_gfa(
            sg,
            sources,
            read_to_unitig,
            (opts.max_tip * 2) as u32,
            options.maximum_hang_length,
            options.minimum_overlap_length,
            options.maximum_hang_rate,
            maximum_length,
        );
    }
}
