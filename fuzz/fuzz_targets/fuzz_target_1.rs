#![no_main]

extern crate libfuzzer_sys;
extern crate arbitrary;
extern crate sevec;

use libfuzzer_sys::fuzz_target;
use sevec::Sevec;

#[derive(Debug, arbitrary::Arbitrary)]
enum Op {
    /// Tests [`Sevec::insert`].
    Insert { index: usize, value: u8 },
    // InsertStaticSlice { index: usize, value: &[u8] },
    /// Tests [`Sevec::insert_slice`].
    InsertSlice { index: usize, value: Vec<u8> },
    /// Tests [`Sevec::remove`].
    Remove { index: usize },
    /// Tests [`Sevec::remove_and_copy_slice_from_end`].
    RemoveAndCopySliceFromEnd { amnt: usize },
    /// Tests [`Sevec::remove_range`].
    RemoveRange { start: usize, end: usize },
    /// Removes all elements after a given index.
    /// This is to cover the code path where we call [`Sevec::remove_range`] with something like
    /// `..n`.
    RemoveRangeEndBound { start: usize },
    /// Tests [`Sevec::push`].
    Push { value: u8 },
    /// Tests [`Sevec::push_slice`].
    PushSlice { value: Vec<u8> },
    // Set { index: usize, value: u8 },
    /// Tests [`Sevec::get`].
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

            Op::RemoveAndCopySliceFromEnd { amnt } => {

                let amnt = amnt % (refr.len() + 2);

                trace.push(format!("RemoveRange {{ amnt: {amnt} }}"));

                if let Some(data) = data.remove_and_copy_slice_from_end(amnt) {
                    assert_eq!(data.len(), amnt);
                    for _i in 0..amnt {
                        refr.pop().unwrap();
                    }
                }
                else {
                    assert!(refr.len() <= amnt);
                }

            },

            Op::RemoveRange { start, end } => {

                // Initializes Values
                let start = start % (refr.len() + 2);
                let end = end % (refr.len() + 2);

                trace.push(format!("RemoveRange {{ start: {start}, end: {end} }}"));

                if let Some(()) = data.remove_range(start..=end) {

                    // Removes the start index the correct amount of times.
                    for i in start..=end {
                        refr.remove(start);
                    }

                }
                else {
                    assert!(
                        start >= refr.len() ||
                        end >= refr.len() ||
                        end < start
                    );
                }

            },

            Op::RemoveRangeEndBound { start } => {

                let start = start % (refr.len() + 2);

                trace.push(format!("RemoveRangeEndBound {{ start: {start} }}"));

                if let Some(()) = data.remove_range(start..) {

                    // Removes the start index the correct amount of times.
                    for i in start..refr.len() {
                        refr.remove(start);
                    }

                }
                else {
                    assert!(
                        start >= refr.len()
                    );
                }


            },

            Op::Push { value } => {
                data.push(value);
                refr.push(value);
            },

            Op::PushSlice { value } => {
                data.push_slice(&value);
                refr.extend(value);
            },

            Op::Get { index } => {
                assert_eq!(data.get(index), refr.get(index));
            },

            // IMPLEND

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
