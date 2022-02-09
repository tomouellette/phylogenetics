use std::collections::HashMap;

// UPGMA (Unweighted Pair Group Method with Arithmetic mean)
// ---------------------------------------------------------
//
// Assumptions for validity:
// - D is a a proper, square distance matrix for all samples for s1...sN
// - A binary tree exists that fulfills the molecular-clock assumption and the distances are given in D
// - D is ultrametric; i.e. for for any triplet {a,b,c} then d_{ab} < d_{ac} = {bc}
pub fn upgma(distance_matrix: Vec<Vec<f64>>, sample_names: Vec<String>) -> HashMap<u64, String> {

    // Update current clusters
    let mut phi = vec![];
    for s in &sample_names {
        phi.push(vec![s.clone()]);
    }

    // Initialize cluster index counter
    let mut ind: u64 = sample_names.len() as u64;    

    // Build a hashmap to track cluster membership and distances
    let mut clustering = HashMap::new();
    for s in 0..sample_names.len() {
        clustering.insert(s as u64, sample_names[s].clone());
    }

    let mut d = distance_matrix.clone();
    while ind < 7 as u64 {

        // Initialize storage vectors
        let mut newick_id = String::new();

        // Find d(i,j) with minimum d(i,j) > 0
        let (mut dij, mut ci, mut cj) = (0.0, 0, 0);
        for i in 0..phi.len() {
            for j in i+1..phi.len() {
                if (d[i][j] < dij) | (i == 0) {
                    dij = d[i][j];
                    println!("dij: {}", dij);
                    ci = i;
                    cj = j;
                }
            }
        }

        // Update newick format 
        // Note: UPGMA leads to identical branch lengths/distances for bifurcating pairs
        let mut newick_vector = vec![
            // ith sequence
            "(".to_string(), phi[cj][0].clone(), ":".to_string(), dij.to_string(), ",".to_string(),
            
            // jth sequence
            phi[ci][0].clone(), ":".to_string(), dij.to_string(), ")".to_string()
        ];        
        newick_id = newick_vector.join("");

        // Replace ci and cj with a single cluster in phi
        phi.remove(ci);
        if ci > cj {
            phi.remove(cj as usize);
        } else {
            phi.remove(cj as usize - 1);
        }
        phi.push(vec![newick_id.clone()]);

        // Add cluster to hash map
        clustering.insert(ind, newick_id.clone());
        
        // Recompute distance matrix for remaining clusters
        d = vec![vec![0.0; phi.len()]; phi.len()];
        for i in 0..phi.len() {
            for j in i+1..phi.len() {
                let mut distance = 0.0;
                d[i][j] = distance;
            }
        }

        ind += 1;
    }

    (clustering)
}

fn main() {
    let distance_matrix = vec![vec![1.0, 4.2, 3.0, 2.0], vec![4.0, 1.0, 2.0, 3.0], vec![3.0, 2.0, 1.0, 3.0], vec![2.0, 3.0, 3.0, 1.0]];
    let names = vec!["a".to_string(), "b".to_string(), "c".to_string(), "d".to_string()];
    let phylo = upgma(distance_matrix, names);
    println!("{:?}", phylo);
}