//! Replay of `/tmp/ed_simi_probe.out`: the no-trace (_simi) edit-distance
//! kernels across 64/128/infi thresholds. Ignored by default.

use rshifiasm::alignment::edit_distance::{
    extension_backward_no_trace, extension_forward_no_trace, global_no_trace,
    semi_global_absent_diagonal_no_trace, EditDistanceResult,
};

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u32 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 33) as u32
    }
}
fn random_sequence(lcg: &mut Lcg, length: usize) -> Vec<u8> {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    (0..length).map(|_| BASES[(lcg.next() & 3) as usize]).collect()
}
fn mutate(lcg: &mut Lcg, p: &[u8], rate: f64) -> Vec<u8> {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let cut = (rate * 1000.0) as u32;
    let mut t = Vec::new();
    for &b in p {
        let r = lcg.next() % 1000;
        if r < cut {
            let k = lcg.next() % 3;
            if k == 0 {
                t.push(BASES[((lcg.next() >> 2) & 3) as usize]);
            } else if k == 1 {
                t.push(BASES[((lcg.next() >> 2) & 3) as usize]);
                t.push(b);
            }
        } else {
            t.push(b);
        }
    }
    t
}

fn fresh() -> EditDistanceResult {
    EditDistanceResult {
        threshold: 0, error: i32::MAX,
        pattern_start: 0, pattern_end: 0, text_start: 0, text_end: 0,
        pattern_length: 0, text_length: 0,
    }
}

#[test]
#[ignore]
fn simi_kernels_match_the_reference() {
    let golden = std::fs::read_to_string("/tmp/ed_simi_probe.out").unwrap();
    let mut lines = golden.lines();
    let mut lcg = Lcg(88172645463325252);
    let thres = [5, 15, 31, 32, 40, 63, 64, 80, 100, 150, 200];
    let mut id = 0;
    let mut mismatches = 0;
    for &thre in &thres {
        for rep in 0..6 {
            let pn = 200 + (lcg.next() % 600) as usize;
            let pattern = random_sequence(&mut lcg, pn);
            let rate = 0.02 + 0.04 * rep as f64;
            let text = mutate(&mut lcg, &pattern, rate);
            let tn = text.len();
            let aux = (lcg.next() % (2 * thre as u32 + 1)) as i32;

            let mut checks: Vec<(&str, EditDistanceResult)> = Vec::new();
            let mut r = fresh();
            global_no_trace(&pattern, pn as i32, &text, tn as i32, thre, &mut r);
            checks.push(("G", r));
            let mut r = fresh();
            extension_forward_no_trace(&pattern, pn as i32, &text, tn as i32, thre, &mut r);
            checks.push(("E0", r));
            let mut r = fresh();
            extension_backward_no_trace(&pattern, pn as i32, &text, tn as i32, thre, &mut r);
            checks.push(("E1", r));
            let mut r = fresh();
            semi_global_absent_diagonal_no_trace(&pattern, pn as i32, &text, tn as i32, thre, aux, &mut r);
            checks.push(("S", r));

            for (fn_name, r) in checks {
                let line = lines.next().expect("golden line");
                let expected = format!(
                    "{} id={} pn={} tn={} thre={} aux={} nword={} err={} ps={} pe={} ts={} te={}",
                    fn_name, id, pn, tn, thre, if fn_name == "S" { aux } else { 0 },
                    0, r.error, r.pattern_start, r.pattern_end, r.text_start, r.text_end
                );
                // nword in golden is ez->nword which the simi kernels may not set; compare rest
                let norm = |s: &str| -> String {
                    s.split_whitespace()
                        .filter(|t| !t.starts_with("nword="))
                        .collect::<Vec<_>>()
                        .join(" ")
                };
                if norm(line) != norm(&expected) {
                    mismatches += 1;
                    if mismatches <= 30 {
                        eprintln!("MISMATCH\n  cxx: {}\n  rs : {}", line, expected);
                    }
                }
            }
            id += 1;
        }
    }
    assert_eq!(mismatches, 0, "{} kernel mismatches", mismatches);
}