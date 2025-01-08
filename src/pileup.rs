use log::{debug, info, warn};
use csv::{ByteRecord, ReaderBuilder};
use core::panic;
use std::io::{Read};
use std::collections::VecDeque;
use std::thread::current;
use anyhow::anyhow;
use anyhow::{Result};
use atoi;
use motif_methylation_state::utils::{
    iupac, 
    modtype, 
    motif,
    modtype::ModType,
    strand::Strand,
};

#[derive(Debug, Clone)]
pub struct PileupRecord {
    pub reference: String,
    pub position: usize,
    pub strand: Strand,
    pub mod_type: ModType,
    pub n_mod: u32,
    pub n_valid_cov: u32,
    pub n_canonical: u32,
    pub n_diff: u32,
}

#[derive(Debug)]
pub struct PileupChunk {
    pub reference: String,
    pub records: Vec<PileupRecord>,
}
/// Reader for pileup chunks
pub struct PileupChunkReader<R: Read> {
    reader: csv::Reader<R>,
    buffer: VecDeque<ByteRecord>,
    min_cov: u32,
    pub eof_reached: bool,
}

impl<R: Read> PileupChunkReader<R> {
    /// Creates a new `PileupChunkReader`
    pub fn new(inner: R, min_cov: u32) -> Self {
        let reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(inner);
        Self {
            reader,
            buffer: VecDeque::new(),
            min_cov,
            eof_reached: false,
        }
    }

    /// Reads the next chunk of records grouped by the same reference
    pub fn next_chunk(&mut self) -> Option<PileupChunk> {
        let mut parsed_records = Vec::new();
        let mut current_reference = None;
        let mut record = ByteRecord::new();
        // Process records from the buffer
        while let Some(record) = self.buffer.pop_front() {
            let reference = std::str::from_utf8(record.get(0).unwrap_or(b""))
                .unwrap_or("")
                .to_string();
    
            match &current_reference {
                Some(current_ref) if reference != *current_ref => {
                    self.buffer.push_front(record);
                    break;
                }
                None => {
                    current_reference = Some(reference.clone());
                }
                _ => {}
            }
    
            if let Ok(parsed_record) = parse_and_validate_pileup_record(&record, self.min_cov) {
                parsed_records.push(parsed_record);
            }
        }
        
        // Load the next batch of records
        while let Ok(has_record) = self.reader.read_byte_record(&mut record) {
            if !has_record {
                self.eof_reached = true;
                break;
            }
            let reference = std::str::from_utf8(record.get(0).unwrap_or(b"")).unwrap();
    
            match current_reference {
                Some(current_ref) if reference != current_ref => {
                    self.buffer.push_front(record);
                    break;
                }
                None => {
                    current_reference = Some(reference.to_string());
                }
                _ => {}
            }
    
            if let Ok(parsed_record) = parse_and_validate_pileup_record(&record, self.min_cov) {
                parsed_records.push(parsed_record);
            }
        }
        if parsed_records.is_empty() {
            return None;
        }
        // Get the reference from the first record
        let reference_out = parsed_records[0].reference.clone();
        Some(PileupChunk {
            reference:  reference_out,
            records: parsed_records,
        })
    }

    pub fn load_n_chunks(&mut self, n: usize) -> Option<Vec<PileupChunk>> {
        let mut chunks = Vec::new();
        for _ in 0..n {
            if let Some(chunk) = self.next_chunk() {
                chunks.push(chunk);
            } else {
                break; // Stop if no more chunks are available
            }
        }
    
        if chunks.is_empty() {
            None
        } else {
            Some(chunks)
        }
    }
}
fn parse_and_validate_pileup_record(
    record: &ByteRecord,
    min_cov: u32,
) -> Result<PileupRecord> {
    let n_valid_cov: u32 = atoi::atoi(record.get(9).unwrap_or(b""))
        .ok_or_else(|| anyhow!("Invalid n_valid_cov value"))?;
    if n_valid_cov < min_cov {
        return Err(anyhow!("Coverage below minimum threshold"));
    }
    let reference = std::str::from_utf8(record.get(0).unwrap_or(b""))?.to_string();
    let position = atoi::atoi(record.get(1).unwrap_or(b""))
        .ok_or_else(|| anyhow!("Could not parse pileup position value"))?;
    let strand = std::str::from_utf8(record.get(5).unwrap_or(b""))?
        .parse::<Strand>()
        .map_err(|_| anyhow!("Could not parse pileup strand value"))?;
    let mod_type = std::str::from_utf8(record.get(3).unwrap_or(b""))?
        .parse::<ModType>()
        .map_err(|_| anyhow!("Could not parse pileup mod_type value"))?;
    let n_mod = atoi::atoi(record.get(11).unwrap_or(b""))
        .ok_or_else(|| anyhow!("Could not parse pileup n_mod value"))?;
    let n_canonical = atoi::atoi(record.get(12).unwrap_or(b""))
        .ok_or_else(|| anyhow!("Could not parse pileup n_canonical value"))?;
    let n_diff = atoi::atoi(record.get(17).unwrap_or(b""))
        .ok_or_else(|| anyhow!("Could not parse pileup n_diff value"))?;
    let pileup_record = PileupRecord {
        reference: reference,
        position: position,
        strand: strand,
        mod_type: mod_type,
        n_mod: n_mod,
        n_valid_cov: n_valid_cov,
        n_canonical: n_canonical,
        n_diff: n_diff,
    };
    Ok(pileup_record)
}



impl PileupRecord {
    fn is_valid(&self, min_cov: u32) -> bool {
        if self.n_valid_cov < min_cov {
            return false;
        }
        true
    }
    fn reference(&self) -> &str {
        &self.reference
    }
}
