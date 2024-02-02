#![allow(unused, non_snake_case)]

#[macro_use]
extern crate lazy_static;

use parking_lot::{MutexGuard, Mutex};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use serde::{Deserialize, Serialize};
use serde_json::from_str;
use std::collections::{HashMap, HashSet};
use std::ops::DerefMut;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::Duration;
use std::{panic, path::PathBuf, thread};
use std::cell::RefCell;

lazy_static! {
    static ref IUPAC_CODES:HashMap<u8, Vec<u8>> = HashMap::from([
            (b'R', vec![b'A', b'G']),
            (b'Y', vec![b'C', b'T']),
            (b'S', vec![b'G', b'C']),
            (b'W', vec![b'A', b'T']),
            (b'K', vec![b'G', b'T']),
            (b'M', vec![b'A', b'C']),
            (b'B', vec![b'C', b'G', b'T']),
            (b'D', vec![b'A', b'G', b'T']),
            (b'H', vec![b'A', b'C', b'T']),
            (b'V', vec![b'A', b'C', b'G']),
            (b'N', vec![b'A', b'C', b'G', b'T'])
        ]);
    static ref VEC_PEPS: Vec<char> = vec![
        'A', 'L', 'W', 'Q', 'Y', 'E', 'C', 'D', 'F', 'G', 'H', 'I', 'M', 'K', 'P', 'R',
        'S', 'V', 'N', 'T', '*', '-', 'B', 'J', 'Z', 'X',
    ];
    // static ref BASES: Vec<char> = vec!['A', 'T', 'G', 'C', 'U', 'N'];

    static ref DICT_TABLE: HashMap<i32, HashMap<char, Vec<String>>> = HashMap::from([
//    let one:  HashMap<char,Vec<String>> = HashMap::from([
        (1, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TAG".to_string(), "TGA".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),

    //let two: HashMap<char, Vec<String>> = HashMap::from([
        (2, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TAG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGA".to_string(), "TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string()]),
            ('M', vec!["ATA".to_string(), "ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
            ])),

        (3, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGA".to_string(), "TGG".to_string()]),
            ('T', vec!["CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string()]),
            ('M', vec!["ATA".to_string(), "ATG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "ATT".to_string(), "ATC".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),

        ])),
        (4, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGA".to_string(), "TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (5, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGA".to_string(), "TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string()]),
            ('M', vec!["ATA".to_string(), "ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (6, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('Q', vec!["TAA".to_string(), "TAG".to_string(), "CAA".to_string(), "CAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('*', vec!["TGA".to_string()]),
            ('W', vec!["TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["TAA".to_string(), "TAG".to_string(), "CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),

        (9, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGA".to_string(), "TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string(), "AAA".to_string()]),
            ('K', vec!["AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "AAA".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (10, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string(), "TGA".to_string()]),
            ('W', vec!["TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (11, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TAG".to_string(), "TGA".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (12, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "CTG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TAG".to_string(), "TGA".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (13, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGA".to_string(), "TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string()]),
            ('M', vec!["ATA".to_string(), "ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('G', vec!["AGA".to_string(), "AGG".to_string(), "GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (14, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string(), "TAA".to_string()]),
            ('*', vec!["TAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGA".to_string(), "TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string(), "AAA".to_string()]),
            ('K', vec!["AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "AAA".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (15, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TGA".to_string()]),
            ('Q', vec!["TAG".to_string(), "CAA".to_string(), "CAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["TAG".to_string(), "CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (16, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "TAG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TGA".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "TAG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (21, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGA".to_string(), "TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string()]),
            ('M', vec!["ATA".to_string(), "ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string(), "AAA".to_string()]),
            ('K', vec!["AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "AAA".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (22, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "TAG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('*', vec!["TCA".to_string(), "TAA".to_string(), "TGA".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "TAG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (23, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('*', vec!["TTA".to_string(), "TAA".to_string(), "TAG".to_string(), "TGA".to_string()]),
            ('L', vec!["TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (24, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string(), "AGA".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGA".to_string(), "TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string(), "AGG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (25, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('G', vec!["TGA".to_string(), "GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('W', vec!["TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (26, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TAG".to_string(), "TGA".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGG".to_string()]),
            ('A', vec!["CTG".to_string(), "GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (27, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('Q', vec!["TAA".to_string(), "TAG".to_string(), "CAA".to_string(), "CAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGA".to_string(), "TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["TAA".to_string(), "TAG".to_string(), "CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (28, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('Q', vec!["TAA".to_string(), "TAG".to_string(), "CAA".to_string(), "CAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGA".to_string(), "TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["TAA".to_string(), "TAG".to_string(), "CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (29, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string(), "TAA".to_string(), "TAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('*', vec!["TGA".to_string()]),
            ('W', vec!["TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (30, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('E', vec!["TAA".to_string(), "TAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('*', vec!["TGA".to_string()]),
            ('W', vec!["TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["TAA".to_string(), "TAG".to_string(), "GAA".to_string(), "GAG".to_string(), "CAA".to_string(), "CAG".to_string()]),
        ])),
        (31, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('E', vec!["TAA".to_string(), "TAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGA".to_string(), "TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["TAA".to_string(), "TAG".to_string(), "GAA".to_string(), "GAG".to_string(), "CAA".to_string(), "CAG".to_string()]),
        ])),
        (32, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string()]),
            ('*', vec!["TAA".to_string(), "TGA".to_string()]),
            ('W', vec!["TAG".to_string(), "TGG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string(), "AGA".to_string(), "AGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),
        (33, HashMap::from([
            ('F', vec!["TTT".to_string(), "TTC".to_string()]),
            ('L', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string()]),
            ('S', vec!["TCT".to_string(), "TCC".to_string(), "TCA".to_string(), "TCG".to_string(), "AGT".to_string(), "AGC".to_string(), "AGA".to_string()]),
            ('Y', vec!["TAT".to_string(), "TAC".to_string(), "TAA".to_string()]),
            ('*', vec!["TAG".to_string()]),
            ('C', vec!["TGT".to_string(), "TGC".to_string()]),
            ('W', vec!["TGA".to_string(), "TGG".to_string()]),
            ('P', vec!["CCT".to_string(), "CCC".to_string(), "CCA".to_string(), "CCG".to_string()]),
            ('H', vec!["CAT".to_string(), "CAC".to_string()]),
            ('Q', vec!["CAA".to_string(), "CAG".to_string()]),
            ('R', vec!["CGT".to_string(), "CGC".to_string(), "CGA".to_string(), "CGG".to_string()]),
            ('I', vec!["ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('M', vec!["ATG".to_string()]),
            ('T', vec!["ACT".to_string(), "ACC".to_string(), "ACA".to_string(), "ACG".to_string()]),
            ('N', vec!["AAT".to_string(), "AAC".to_string()]),
            ('K', vec!["AAA".to_string(), "AAG".to_string(), "AGG".to_string()]),
            ('V', vec!["GTT".to_string(), "GTC".to_string(), "GTA".to_string(), "GTG".to_string()]),
            ('A', vec!["GCT".to_string(), "GCC".to_string(), "GCA".to_string(), "GCG".to_string()]),
            ('D', vec!["GAT".to_string(), "GAC".to_string()]),
            ('E', vec!["GAA".to_string(), "GAG".to_string()]),
            ('G', vec!["GGT".to_string(), "GGC".to_string(), "GGA".to_string(), "GGG".to_string()]),
            ('X', vec![]),
            ('B', vec!["AAT".to_string(), "AAC".to_string(), "GAT".to_string(), "GAC".to_string()]),
            ('J', vec!["TTA".to_string(), "TTG".to_string(), "CTT".to_string(), "CTC".to_string(), "CTA".to_string(), "CTG".to_string(), "ATT".to_string(), "ATC".to_string(), "ATA".to_string()]),
            ('Z', vec!["CAA".to_string(), "CAG".to_string(), "GAA".to_string(), "GAG".to_string()]),
        ])),

                ]);

}
fn recurse(triplet: &[u8], working: &mut [u8], index: usize, output: &mut HashSet<String>) {
    if index == triplet.len() {
        // println!("{}", String::from_utf8(working.to_vec()).unwrap());
        output.insert(String::from_utf8(working.to_vec()).unwrap());
        return;
    }
    let this_char = triplet[index];
    recurse(triplet, working, index + 1, output);
    if IUPAC_CODES.contains_key(&this_char) {
        for character in IUPAC_CODES.get(&this_char).unwrap() {
            working[index] = *character;
            recurse(triplet, working, index+1, output);
        }
    }
    return;
}

pub fn make_iupac_set(triplet: &[u8]) -> HashSet<String> {
    let mut output = HashSet::new();
    let mut working = [0_u8;3];
    for (i, &byte) in triplet.iter().enumerate(){
        working[i] = byte;
    }
    recurse(triplet, &mut working, 0_usize, &mut output);
    output
}

#[pyfunction]
pub fn attempt_iupac_substitution(original_triplet: &str, taxa: Vec<String>) -> Option<String> {
    let original_bytes = original_triplet.as_bytes();
    let possible_subs = make_iupac_set(original_bytes);
    // println!("searching for {}", original_triplet);
    for made in &possible_subs {
        // println!("{}", made);
    }

    for triplet in &taxa {
        if possible_subs.contains(triplet) {
            return Some(triplet.to_string());
        }
    }
    None
}

#[derive(Clone)]
struct AminoAcidTranslator(
    (String, String), 
    (String, String),
    (RefCell<bool>, String),
);

impl AminoAcidTranslator {
    pub fn do_checks(&self) {
        let AminoAcidTranslator((aa_header, aa), (nt_header, nt), _) = self;

        if aa_header != nt_header {
            self.error_out(format!(
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

            self.error_out(format!(
                "{} is larger than {} {}",
                longer_shorter.0, longer_shorter.1, diff
            ));
        }
    }

    pub fn streamline(&mut self) {
        let AminoAcidTranslator((header, amino_acid), (_, nucleotide), _) = self;

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

    fn error_out(&self, message: String) {
        let AminoAcidTranslator((header, _), _, (dont_skip, file_stem)) = self;      

        match dont_skip.clone().into_inner() {
            true => {
                let mut dont_skip_deref = dont_skip.borrow_mut();

                *dont_skip_deref = false;

                println!(
                    "\n===ERROR CAUGHT IN FILE {} AND HEADER {}:\n {}\n===",
                    file_stem, header, message
                );
            },
            false => (),
        }

    }

    fn error_out_mismatch(&self, aa_counter: usize) {
        let AminoAcidTranslator((header, amino_acid), (_, compare_dna), _) = self;
        let start = std::cmp::max(10, aa_counter) - 10;

        let stop = std::cmp::min(amino_acid.len() - 10, aa_counter) + 10;
        let percent = aa_counter as f64/amino_acid.len() as f64;
        let compare_index = (percent * compare_dna.len() as f64) as usize;

        let compare_start = std::cmp::max(30, compare_index) - 30;
        let compare_end = std::cmp::min(compare_dna.len() - 30, compare_index) + 30 ;

        self.error_out(format!(
            r#" 
                ======
                MISMATCH ERROR:
                The following Amino Acid failed to match with its source Nucleotide pair at aa site {}.
                Amino Acid: `{}`,
                ======
                Source Nucleotide: `{}`,
                =======
            "#,
            aa_counter, &amino_acid[start..stop], &compare_dna[compare_start..compare_end]
        ));
    }

    pub fn reverse_translate_and_compare(&self, gene_table: &HashMap<char, Vec<String>>) -> String {
        let AminoAcidTranslator((header, amino_acid), (_, compare_dna), _) = self;

        let mut compare_triplets = (0..compare_dna.len())
            .step_by(3)
            .map(|i| compare_dna[i..i + 3].to_string())
            .into_iter();
        let mut aa_counter: usize = 0;
        amino_acid
            .chars()
            .enumerate()
            .map(|(aa_index, aa)| {
                // aa_counter += 1;
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
                                                        // println!{"match failed, attemptint to rescue {}", &original_triplet}
                                                        match attempt_iupac_substitution(&original_triplet, taxa.clone()) {
                                                            Some(t) => return t,
                                                            None => { self.error_out_mismatch(aa_index);
                                                                return "".to_string();
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }                                        
                                    }
                                    None => {
                                        self.error_out(
                                            "Genetic table does not have the pep. Perhaps you've chosen the wrong table index?".to_string()
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
// There is something weird going on with the second seqs tuple.
// It used to run with two elements, but now I get a value error because its taking 3,
// but the tuple hasnt changed. Should I make two versions of the function?
#[pyfunction]
pub fn pn2codon(
    file_steem: String,
    // gene_table: HashMap<char, Vec<String>>,
    table_num: i32,
    seqs: HashMap<String, ((String, String), (i32, String, String))>
) -> String {    
    let mut dont_skip = RefCell::new(true);
    let gene_table = DICT_TABLE.get(&table_num).unwrap();
    let file = String::from_iter(
        seqs
            .iter()
            .take_while(|_| dont_skip.clone().into_inner())
            .map(|(header, ((aa_header, aa), (_ ,nt_header, nt)))| {
                if !dont_skip.clone().into_inner() {
                    println!("{}", dont_skip.clone().into_inner());
                }
                // let gene_table = DICT_TABLE.get(table_key).unwrap();
                let mut amino_acid = AminoAcidTranslator(
                    (aa_header.clone(), aa.clone()),
                    (nt_header.clone(), nt.clone()),
                    (dont_skip.clone(), file_steem.clone())
                );
                amino_acid.streamline();
                amino_acid.do_checks();

                let mut codon = amino_acid.reverse_translate_and_compare(&gene_table);
                let mut h_clone = header.clone();

                codon.push('\n');
                h_clone.push('\n');

                vec![h_clone, codon]
            })
            .flatten()
    );

    file
}

#[pyfunction]
pub fn pn2codon_original_args(
    file_steem: String,
    gene_table: HashMap<char, Vec<String>>,
    seqs: HashMap<String, ((String, String), (String, String))>
) -> String {
    let mut dont_skip = RefCell::new(true);

    let file = String::from_iter(
        seqs
            .iter()
            .take_while(|_| dont_skip.clone().into_inner())
            .map(|(header, ((aa_header, aa), (nt_header, nt)))| {
                if !dont_skip.clone().into_inner() {
                    println!("{}", dont_skip.clone().into_inner());
                }

                let mut amino_acid = AminoAcidTranslator(
                    (aa_header.clone(), aa.clone()),
                    (nt_header.clone(), nt.clone()),
                    (dont_skip.clone(), file_steem.clone())
                );
                amino_acid.streamline();
                amino_acid.do_checks();

                let mut codon = amino_acid.reverse_translate_and_compare(&gene_table);
                let mut h_clone = header.clone();

                codon.push('\n');
                h_clone.push('\n');

                vec![h_clone, codon]
            })
            .flatten()
    );

    file
}

#[pymodule]
fn pr2codon(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(pn2codon, m)?)?;
    m.add_function(wrap_pyfunction!(pn2codon_original_args, m)?)?;
    m.add_function(wrap_pyfunction!(attempt_iupac_substitution, m)?)?;
    Ok(())
}
