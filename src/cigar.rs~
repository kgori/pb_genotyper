use rust_htslib::bam::Record;

pub(crate) fn cigar_end_pos(record: &Record) -> i64 {
    let mut pos = record.pos(); // 0-based position of the read

    for cigar in record.cigar().iter() {
        let len = cigar.len() as i64;
        match cigar.char() {
            'M' | '=' | 'X' | 'D' | 'N' => pos += len, // Consider 'M', '=', 'X', 'D', and 'N' as alignment matches
            'I' | 'S' | 'H' | 'P' => {} // Insertion to the reference, soft clipping, hard clipping, and padding (ignored for alignment)
            _ => {}
        }
    }
    pos
}

// Function to check if there's an overlap between the read (based on CIGAR) and an interval
// The interval is 0-based, half-open, i.e. [start, end)
pub(crate) fn check_cigar_overlap(record: &Record, interval_start: i64, interval_end: i64) -> bool {
    let mut pos = record.pos(); // 0-based position of the read

    for cigar in record.cigar().iter() {
        let len = cigar.len() as i64;
        match cigar.char() {
            'M' | '=' | 'X' => {
                // Consider 'M', '=', and 'X' as alignment matches
                if pos < interval_end && pos + len > interval_start {
                    // There is an overlap
                    return true;
                }
                pos += len; // Move past this segment for the next CIGAR element
            }
            'D' | 'N' => pos += len, // Deletion or skipped region from the reference
            'I' | 'S' | 'H' | 'P' => {} // Insertion to the reference, soft clipping, hard clipping, and padding (ignored for alignment)
            _ => {}
        }
    }
    false // No overlap found
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::{
        bam,
        bam::record::{Cigar, CigarString},
    };

    // Helper function to create a mock BAM record with a given sequence of CIGAR operations and starting position
    fn mock_record(cigar_ops: Vec<(char, u32)>, start_pos: i64) -> bam::Record {
        let mut record = bam::Record::new();
        let cigar = CigarString(
            cigar_ops
                .iter()
                .map(|&(op, len)| match op {
                    'M' => Cigar::Match(len),
                    'I' => Cigar::Ins(len),
                    'D' => Cigar::Del(len),
                    'S' => Cigar::SoftClip(len),
                    'H' => Cigar::HardClip(len),
                    'N' => Cigar::RefSkip(len),
                    '=' => Cigar::Equal(len),
                    'X' => Cigar::Diff(len),
                    _ => panic!("Unsupported CIGAR operation"),
                })
                .collect(),
        );
        // Set dummy sequence and qualities to allow setting CIGAR
        let seq_len = cigar
            .iter()
            .filter_map(|c| match c {
                Cigar::Match(l)
                | Cigar::Equal(l)
                | Cigar::Diff(l)
                | Cigar::SoftClip(l)
                | Cigar::Ins(l) => Some(*l as usize),
                _ => None,
            })
            .sum();
        let seq = vec![0; seq_len]; // Dummy sequence
        let qual = vec![255; seq_len]; // Dummy quality scores
        let qname = vec![b'A' as u8]; // Dummy read name
        record.set(&qname, Some(&cigar), &seq, &qual);
        record.set_pos(start_pos);

        // Since there's no direct set_cigar, we manipulate data directly for testing purpose
        // This is not typical usage, and should be used with caution outside of testing
        record.cache_cigar();
        record
    }

    #[test]
    fn test_cigar_end_pos_simple() {
        let record = mock_record(vec![('M', 50)], 100); // 50 matches from position 100
        assert_eq!(cigar_end_pos(&record), 150); // End position of the match
    }

    #[test]
    fn test_cigar_end_pos_refskip() {
        let record = mock_record(vec![('M', 10), ('N', 50), ('M', 10)], 100); // Complex CIGAR starting at position 20
        assert_eq!(cigar_end_pos(&record), 170); // End position of the match
    }

    #[test]
    fn test_cigar_end_pos_rnaseq() {
        let record = mock_record(vec![('M', 85), ('N', 24_899), ('M', 16)], 61_820_205);
        assert_eq!(cigar_end_pos(&record), 61_845_205); // End position of the match
    }

    #[test]
    fn test_overlap_full_match() {
        let record = mock_record(vec![('M', 50)], 100); // 50 matches from position 100
        assert!(check_cigar_overlap(&record, 120, 130)); // Overlaps within the match
    }

    #[test]
    fn test_overlap_boundaries_single_match() {
        let record = mock_record(vec![('M', 1)], 100); // 50 matches from position 100
        assert!(!check_cigar_overlap(&record, 90, 100)); // Closest non-overlapping position
        assert!(check_cigar_overlap(&record, 90, 101)); // Overlaps at the start of the match
        assert!(check_cigar_overlap(&record, 100, 110)); // Overlaps at the end of the match
        assert!(!check_cigar_overlap(&record, 101, 110)); // Closest non-overlapping position
    }

    #[test]
    fn test_overlap_boundaries_region_match() {
        let record = mock_record(vec![('M', 10)], 100); // 50 matches from position 100
        assert!(!check_cigar_overlap(&record, 90, 100)); // Closest non-overlapping position
        assert!(check_cigar_overlap(&record, 90, 101)); // Overlaps at the start of the match
        assert!(check_cigar_overlap(&record, 109, 120)); // Overlaps at the end of the match
        assert!(!check_cigar_overlap(&record, 110, 120)); // Closest non-overlapping position
    }

    #[test]
    fn test_overlap_no_match() {
        let record = mock_record(vec![('M', 50)], 100); // 50 matches from position 100
        assert!(!check_cigar_overlap(&record, 151, 160)); // No overlap, outside the match
    }

    #[test]
    fn test_overlap_partial_match() {
        let record = mock_record(vec![('M', 30), ('S', 20)], 100); // 30 matches then 20 soft clipped from position 100
        assert!(check_cigar_overlap(&record, 120, 130)); // Overlaps within the match
        assert!(!check_cigar_overlap(&record, 130, 140)); // No overlap, within soft clip
    }

    #[test]
    fn test_overlap_with_insertion() {
        let record = mock_record(vec![('M', 20), ('I', 10), ('M', 20)], 100); // 20 matches, 10 insertion, then 20 matches
        assert!(check_cigar_overlap(&record, 115, 125)); // Overlaps across the insertion
    }

    #[test]
    fn test_overlap_with_deletion() {
        let record = mock_record(vec![('M', 20), ('D', 10), ('M', 30)], 100); // 20 matches, 10 deletion, then 30 matches
        assert!(check_cigar_overlap(&record, 115, 125)); // Overlaps before the deletion
        assert!(check_cigar_overlap(&record, 130, 140)); // Overlaps after the deletion
    }

    #[test]
    fn test_overlap_with_match_at_start() {
        let record = mock_record(vec![('M', 100)], 0); // Matches from position 0 to 100
        assert!(check_cigar_overlap(&record, 1, 50)); // Interval starts at start of read
    }

    #[test]
    fn test_overlap_with_match_at_end() {
        let record = mock_record(vec![('M', 100)], 50); // Matches from position 50 to 150
        assert!(check_cigar_overlap(&record, 100, 150)); // Interval overlaps end of read
    }

    #[test]
    fn test_no_overlap_before_read() {
        let record = mock_record(vec![('M', 100)], 100); // Matches from position 100 to 200
        assert!(!check_cigar_overlap(&record, 1, 50)); // Interval is completely before read
    }

    #[test]
    fn test_no_overlap_after_read() {
        let record = mock_record(vec![('M', 100)], 0); // Matches from position 0 to 100
        assert!(!check_cigar_overlap(&record, 150, 200)); // Interval is completely after read
    }

    #[test]
    fn test_overlap_with_complex_cigar() {
        let record = mock_record(vec![('M', 20), ('I', 10), ('D', 10), ('M', 60)], 20); // Complex CIGAR starting at position 20
        assert!(check_cigar_overlap(&record, 25, 45)); // Overlaps the first 'M' and 'D'
        assert!(check_cigar_overlap(&record, 45, 55)); // Overlaps the 'D' and second 'M'
    }

    #[test]
    fn test_no_overlap_with_complex_cigar() {
        let record = mock_record(vec![('M', 20), ('I', 10), ('D', 10), ('M', 60)], 20); // Complex CIGAR starting at position 20
        assert!(!check_cigar_overlap(&record, 5, 19)); // Before the read starts
                                                       // assert!(!check_cigar_overlap(&record, 95, 120)); // After the read ends, considering 'D'
        assert!(!check_cigar_overlap(&record, 40, 50)); // Deleted part of read
    }

    #[test]
    fn test_overlap_across_entire_read() {
        let record = mock_record(vec![('M', 100)], 50); // Entire read is a match from position 50
        assert!(check_cigar_overlap(&record, 1, 200)); // Interval covers entire read
    }

    #[test]
    fn test_overlap_single_match_between_refskips() {
        let record = mock_record(
            vec![('M', 10), ('N', 50), ('M', 1), ('N', 50), ('M', 10)],
            100,
        ); // 50 matches between two skipped regions
        assert!(check_cigar_overlap(&record, 150, 170)); // Overlaps within the match
    }

    #[test]
    fn test_overlap_rnaseq_split_read() {
        let record = mock_record(
            vec![('S', 4), ('M', 45), ('N', 25_995), ('M', 26)],
            61_339_734,
        );
        assert!(check_cigar_overlap(&record, 61_339_765, 61_339_999)); // Overlaps within the match
        assert!(check_cigar_overlap(&record, 61_365_782, 61_365_999)); // Overlaps within the match
    }

    #[test]
    fn test_no_overlap_rnaseq_split_read() {
        let record = mock_record(
            vec![('S', 4), ('M', 45), ('N', 25_995), ('M', 26)],
            61_339_734,
        );
        assert!(!check_cigar_overlap(&record, 61_344_619, 61_344_735)); // Overlaps within the match
    }

    #[test]
    fn test_overlap_rnaseq_read_pair() {
        let first_in_pair = mock_record(
            vec![('S', 1), ('M', 76), ('N', 3329), ('M', 24)],
            61_826_914,
        );
        let second_in_pair = mock_record(vec![('S', 2), ('M', 99)], 61_845_189);
        assert!(
            check_cigar_overlap(&first_in_pair, 61_845_189, 61_845_296)
                || check_cigar_overlap(&first_in_pair, 61_830_319, 61_830_481)
        ); // Overlaps within the match
        assert!(
            check_cigar_overlap(&second_in_pair, 61_845_189, 61_845_296)
                || check_cigar_overlap(&second_in_pair, 61_830_319, 61_830_481)
        ); // Overlaps within the match
    }
}
