use wireframes_core::Mesh;

/// Combinatorial closure: every directed edge once with its reverse, and one
/// fan per vertex. Geometric checks are reported separately (rounding to f64
/// can create slivers the exact arrangement did not have).
pub fn closed_manifold(m: &Mesh) -> bool {
    use std::collections::HashMap;
    let mut dir: HashMap<(u32, u32), u32> = HashMap::new();
    for t in m.triangles() {
        if t[0] == t[1] || t[1] == t[2] || t[0] == t[2] {
            return false;
        }
        for k in 0..3 {
            *dir.entry((t[k], t[(k + 1) % 3])).or_default() += 1;
        }
    }
    if dir
        .iter()
        .any(|(&(a, b), &n)| n != 1 || dir.get(&(b, a)) != Some(&1))
    {
        return false;
    }
    let mut link: HashMap<u32, HashMap<u32, u32>> = HashMap::new();
    for t in m.triangles() {
        for k in 0..3 {
            link.entry(t[k])
                .or_default()
                .insert(t[(k + 1) % 3], t[(k + 2) % 3]);
        }
    }
    link.values().all(|next| {
        let start = *next.keys().min().unwrap();
        let (mut c, mut n) = (next[&start], 1);
        while c != start && n <= next.len() {
            c = next[&c];
            n += 1;
        }
        n == next.len()
    })
}
