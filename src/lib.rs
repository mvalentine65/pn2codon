#![allow(unused, non_snake_case)]

#[macro_use]
extern crate lazy_static;

use parking_lot::{MutexGuard, Mutex};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use serde::{Deserialize, Serialize};
use serde_json::from_str;
use std::collections::HashMap;
use std::ops::DerefMut;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::Duration;
use std::{panic, path::PathBuf, thread};

static DO_SKIP: AtomicBool = AtomicBool::new(false);


lazy_static! {
    static ref FILE_STEM: Mutex<String> = Mutex::new("".to_string());
    static ref ERROR_HEADER: Mutex<String> = Mutex::new("".to_string());
    static ref VEC_PEPS: Vec<char> = vec![
        'A', 'L', 'W', 'Q', 'Y', 'E', 'C', 'D', 'F', 'G', 'H', 'I', 'M', 'K', 'P', 'R',
        'S', 'V', 'N', 'T', '*', '-', 'B', 'J', 'Z', 'X',
    ];
    static ref BASES: Vec<char> = vec!['A', 'T', 'G', 'C', 'U', 'N'];
}

fn error_out_program(message: &str) {
    println!("Function must error out prematurely:");
    DO_SKIP.store(true, Ordering::Relaxed);
    let file_stem = FILE_STEM.lock();
    let error_header = ERROR_HEADER.lock();


    println!(
        "\n===ERROR CAUGHT IN FILE {} AND HEADER {}: {}\n\n Sleeping for 10 milliseconds to\n===",
        file_stem, error_header, message
    );

    thread::sleep(Duration::from_millis(10));
}

struct AminoAcidTranslator((String, String), (String, String));

impl AminoAcidTranslator {
    pub fn do_checks(&mut self) {
        let AminoAcidTranslator((aa_header, aa), (nt_header, nt)) = self;

        let mut htmx = ERROR_HEADER.lock();
        *htmx = aa_header.clone();
        MutexGuard::unlock_fair(htmx);

        if aa_header != nt_header {
            error_out_program(&format!(
                "AA header -> {} is not the same as NT header -> {}",
                aa_header, nt_header
            ));
        }

        let len_aa = aa.len();
        let len_nt = nt.len();
        let aa_filt_mul = aa.chars().filter(|c| *c != '-').count() * 3;

        if len_nt != aa_filt_mul {
            let longer_shorter = match aa_filt_mul > len_nt {
                true => (
                    format!("(AA -> {})", aa_header),
                    format!("(NT -> {})", nt_header),
                ),
                false => (
                    format!("(NT -> {})", nt_header),
                    format!("(AA -> {})", aa_header),
                ),
            };

            let diff = {
                let num_marker = match aa_filt_mul > len_nt {
                    true => ((aa_filt_mul - len_nt) / 3, "PEP char(s)"),
                    false => ((len_nt - aa_filt_mul) / 3, "NT triplet(s)"),
                };
                format!("with a difference of {} {}", num_marker.0, num_marker.1)
            };

            error_out_program(&format!(
                "{} is larger than {} {}",
                longer_shorter.0, longer_shorter.1, diff
            ));
        }
    }

    pub fn streamline(&mut self) {
        let AminoAcidTranslator((header, amino_acid), (_, nucleotide)) = self;

        let mut amino_acid_trimmed = amino_acid.trim().to_uppercase();

        let mut amino_acid_filtered = String::new();

        amino_acid_trimmed.char_indices().for_each(|(i, c)| {
            match !VEC_PEPS.contains(&c)
            {
                true => {
                    amino_acid_filtered.push('X');
                }
                false => amino_acid_filtered.push(c),
            }
        });

        *amino_acid = amino_acid_filtered;

        *nucleotide = nucleotide.replace("-", "").replace(".", "");
    }

    fn error_out(&self) {
        let AminoAcidTranslator((header, amino_acid), (_, compare_dna)) = self;

        error_out_program(&format!(
            r#" 
                ======
                MISMATCH ERROR:
                The following Amino Acid failed to match with its source Nucleotide pair.

                Amino Acid: `{}`,
                ======
                Source Nucleotide: `{}`,
                =======
            "#,
            amino_acid, compare_dna
        ));
    }

    pub fn reverse_translate_and_compare(&self, gene_table: &HashMap<char, Vec<String>>) -> String {
        let AminoAcidTranslator((header, amino_acid), (_, compare_dna)) = self;

        let mut compare_triplets = (0..compare_dna.len())
            .step_by(3)
            .map(|i| compare_dna[i..i + 3].to_string())
            .into_iter();

        amino_acid
            .chars()
            .map(|aa| {
                match aa == '-' {
                    true => {
                        return "---".to_string() 
                    },                        
                    false => {
                        match aa.is_ascii_digit() {
                            true => {
                                let d = aa.to_digit(110).unwrap();

                                return ".".repeat(d as usize).to_string()
                            }
                            false => {
                                let mut taxa_triplets = gene_table.get(&aa);

                                match taxa_triplets {
                                    Some(taxa) => {
                                        let mut taxa_mut = taxa.clone();
                                                                           
                                        let original_triplet = compare_triplets.next().unwrap();

                                        match original_triplet.contains('N') || aa == 'X' {
                                            true => { 
                                                return original_triplet;                                                    
                                            },
                                            false => {
                                                taxa_mut.retain(|s| s == &original_triplet);

                                                match taxa_mut.get(0) {
                                                    Some(t) => {
                                                        return t.clone()
                                                    },
                                                    None => {
                                                        self.error_out();
                                                        return "".to_string();
                                                    }
                                                }
                                            }
                                        }                                        
                                    }
                                    None => {
                                        error_out_program(
                                            "Genetic table does not have the pep. Perhaps you've chosen the wrong table index?"
                                        );
                                        return "".to_string();
                                    }
                                }
                            }
                        }
                    }
                }
            })
            .collect::<Vec<String>>()
            .join("")
    }
}

#[pyfunction]
pub fn pn2codon(
    file_steem: String,
    gene_table: HashMap<char, Vec<String>>,
    amino_seqs: Vec<(String, String)>,
    nuc_seqs: Vec<(String, String)>,
) -> String {
    let mut fstem_mutex = FILE_STEM.lock();
    *fstem_mutex = file_steem;
    MutexGuard::unlock_fair(fstem_mutex);

    let aa_seq_len = amino_seqs.len();
    let nt_seq_len = nuc_seqs.len();

    if aa_seq_len != nt_seq_len {
        let longer_shorter = match aa_seq_len > nt_seq_len {
            true => ("AA", "NT"),
            false => ("NT", "AA"),
        };

        let diff = match aa_seq_len > nt_seq_len {
            true => ((aa_seq_len as isize) - (aa_seq_len as isize)).abs(),
            false => ((nt_seq_len as isize) - (aa_seq_len as isize)).abs(),
        };

        println!(
            "Length of the {} sequence is longer than the length of {} sequence by a number of {}.",
            longer_shorter.0, longer_shorter.1, diff
        );

        return "".to_string();
    }
   
    let file = String::from_iter(amino_seqs
        .iter()
        .cloned()
        .zip(nuc_seqs.iter().cloned())
        .take_while(|_| !DO_SKIP.load(Ordering::Relaxed))
        .map(|((aa_header, aa), (nt_header, nt))| {
            let mut amino_acid = AminoAcidTranslator(
                (aa_header.clone(), aa.clone()),
                (nt_header.clone(), nt.clone()),
            );
            amino_acid.do_checks();
            amino_acid.streamline();

            let mut codon = amino_acid.reverse_translate_and_compare(&gene_table);
            let mut h_clone = aa_header.clone();

            codon.push('\n');
            h_clone.push('\n');

            vec![h_clone, codon]
        })
        .flatten()
    );

    DO_SKIP.store(false, Ordering::Relaxed);

    file   
}

#[pymodule]
fn pro2codon(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pn2codon, m)?)?;
    Ok(())
}
