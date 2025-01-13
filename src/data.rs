use ahash::{HashMap, HashMapExt};
use crate::pileup::PileupChunk;
use crate::sequence::Contig;







pub struct GenomeWorkSpaceBuilder {
    pub contigs: HashMap<String, Contig>,
}

impl GenomeWorkSpaceBuilder {
    pub fn new() -> Self {
        Self {
            contigs: HashMap::new(),
        }
    }

    pub fn add_contig(&mut self, contig: &Contig) {
        self.contigs.insert(contig.reference.clone(), contig.clone());
    }

    pub fn push_records(&mut self, records: PileupChunk) {
        let reference = records.reference.clone();
        let contig = self
            .contigs
            .get_mut(&reference)
            .expect("Could not find contig");
        contig.add_records(records);
    }

    pub fn build(self) -> GenomeWorkspace {
        GenomeWorkspace {
            contigs: self.contigs,
        }
    }
}


pub struct GenomeWorkspace {
    pub contigs: HashMap<String, Contig>,
}
