#![no_main]

extern crate libfuzzer_sys;
extern crate arbitrary;
extern crate sevec;

use std::panic::catch_unwind;

use libfuzzer_sys::fuzz_target;
use sevec::Sevec;

#[derive(Debug, arbitrary::Arbitrary)]
enum Op {
    Insert { index: usize, value: u8 },
    // InsertStaticSlice { index: usize, value: &[u8] },
    /// Tests [`Sevec::insert_slice`].
    InsertSlice { index: usize, value: Vec<u8> },
    /// Tests [`Sevec::remove`].
    Remove { index: usize },
    RemoveAndCopySliceFromEnd { amnt: usize },
    RemoveRange { start: usize, end: usize },
    /// Removes all elements after a given index.
    /// This is to cover the code path where we call [`Sevec::remove_range`] with something like
    /// `..n`.
    RemoveRangeEndBound { start: usize },
    /// Tests [`Sevec::push`].
    Push { value: u8 },
    PushSlice { value: Vec<u8> },
    Set { index: usize, value: u8 },
    Get { index: usize },
}

fuzz_target!(|ops: Vec<Op>| {

    let mut data: Sevec<u8> = Sevec::new();
    let mut refr: Vec<u8> = Vec::new();

    let mut trace: Vec<String> = Vec::new();

    for op in ops {
        match op {

            Op::Insert { index, value } => {
                let index = index % (refr.len() + 2);
                trace.push(format!("Insert {{ index: {index}, value: {value} }}"));
                if let Some(()) = data.insert_slice(index, &[value]) {
                    refr.insert(index, value);
                }
                else {
                    assert!(index > refr.len());
                }
            },

            Op::InsertSlice { index, value } => {
                let index = index % (refr.len() + 1);
                trace.push(format!("InsertSlice {{ index: {index}, value: {value:?} }}"));
                data.insert_slice(index, &value);
                for i in index..(value.len() + index) {
                    refr.insert(i, value[i - index]);
                }
            },

            Op::Remove { index } => {
                let index = index % (refr.len() + 1);
                trace.push(format!("Remove {{ index: {index} }}"));
                if let Some(()) = data.remove(index) {
                    refr.remove(index);
                }
                else {
                    assert!(index >= refr.len());
                }
            },

            _ => (),

        }

        let mut i = 0;
        for &seg in data.get_refs() {
            for sub_seg in seg {
                if refr.get(i).is_none() || sub_seg != refr.get(i).unwrap() {
                    let sevec_data: Vec<_> = data.into();
                    let trace_str = trace
                        .iter()
                        .map(|v| format!("\n\t{v}"))
                        .collect::<String>()
                        ;

                    panic!("Data Differs!! Sevec: {:?}, Reference: {:?}.\nTrace:{}", sevec_data, refr, trace_str);
                }
                i += 1;
            }
        }

    }


});
