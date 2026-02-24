use std::{ops::{Bound, RangeBounds}, pin::Pin, ptr, sync::Arc};

#[derive(Clone)]
pub struct Sevec<T> {
    /// The underlying data.
    data: Vec<Pin<Arc<[T]>>>,
    /// The ordered references to read through.
    refs: Vec<*const [T]>,
}

impl <T> Sevec<T> {

    /// Creates a new [`Sevec`].
    pub fn new() -> Self {
        return Self {
            data: Vec::new(),
            refs: Vec::new(),
        };
    }

    // Adds a new slice.
    pub fn push_arc_slice(&mut self, value: Pin<Arc<[T]>>) -> () {
        // Gets the reference
        let data_inner_ref = ptr::slice_from_raw_parts(value.as_ptr(), value.len());
        // Adds the data.
        self.data.push(value);
        // Adds the reference.
        self.refs.push(data_inner_ref);
        return ();
    }

    /// Calculates the estimated size of the sevec.
    pub fn size_estimation(&self) -> usize {

        // The size of the inner [`Arc`] type.
        // The [`Pin`] data type doesn't add size but is here because it better reflects the size
        // of the data.
        const ARC_SIZE: usize = size_of::<Pin<Arc<()>>>();

        // this should be equal to sizeof::<usize>() * 2 (len + addr).
        const SLICE_VEC_SIZE: usize = size_of::<Vec<*const [u8]>>();

        let size =
            size_of::<Self>() +
            ARC_SIZE * self.data.len() +
            SLICE_VEC_SIZE * self.refs.len() +
            0
            ;

        return size;

    }

    /// Gets a specified value from both a chunk index and a chunk sub index.
    pub fn get_from_chunk_and_idx(&self, chunk: usize, idx: usize) -> Option<&T> {
        let chunk_slice = self.get_chunk(chunk)?;
        return chunk_slice.get(idx);
    }

    /// Gets a specified chunk as a slice.
    pub fn get_chunk(&self, chunk: usize) -> Option<&[T]> {
        let chunk_ptr = self.refs.get(chunk)?;
        return unsafe { chunk_ptr.as_ref() };
    }

    /// Gets the chunk index of a specified input index.
    /// The first value in the result is the chunk index.
    /// The second value is the sum of all the previous lengths up until that point.
    fn get_chunk_idx(&self, idx: usize) -> Option<(usize, usize)> {

        // Initializes
        let mut total_len = 0;

        // Goes through the references.
        for (i, ref_ptr) in self.refs.iter().enumerate() {

            // Calculates the new length
            let cur_length = ref_ptr.len();
            total_len += cur_length;

            // Checks if we passed it.
            if total_len > idx {
                // Returns the index of the chunk and the sum of previous lengths.
                total_len -= cur_length; // Goes to the start of the selected chunk
                return Some((i, total_len));
            }

        }

        return None;

    }

    /// Gets a reference to some data.
    pub fn get(&self, idx: usize) -> Option<&T> {

        let (chunk_idx, total_len) = self.get_chunk_idx(idx)?;
        let chunk = self.get_chunk(chunk_idx)?;

        // Gets the result
        let final_idx = idx - total_len;
        let res = chunk.get(final_idx)?;
        return Some(res);

    }

}

impl <T: Unpin> Sevec<T> {

    /// Adds a value to the end of the array.
    pub fn push(&mut self, value: T) -> () {

        // Creates new ptr
        let mut arc_ptr = Arc::<[T]>::new_uninit_slice(1);
        // Writes value to it
        let arc_ptr_mut = Arc::get_mut(&mut arc_ptr).unwrap();
        arc_ptr_mut[0].write(value);
        // Removes the "maybe uninit" status
        // This is safe because we literally just initialized the uninitialized value (to the
        // passed value)
        let arc_ptr = unsafe { arc_ptr.assume_init() };
        let arc_ptr = Pin::new(arc_ptr);

        // Pushes the value
        // In theory, we could just hard-code adding the slice length of 1 but I don't really
        // think this would make much of a difference.
        self.push_arc_slice(arc_ptr);
        return ();

    }


    /// Removes the element at a given index.
    pub fn remove(&mut self, idx: usize) -> Option<()> {

        return self.remove_range(idx..=idx);

        /*
        let (chunk_idx, total_len) = self.get_chunk_idx(idx)?;

        // We need to handle several cases.
        // 1. Where it is the only data in the segment.
        // 2. Where it is the start of segmented data.
        // 3. Where it is the end of segmented data.
        // 4. Where it is the middle of segmented data.
        //
        // We are going to start with the simplest case (being case 1).

        let chunk = self.get_chunk(chunk_idx)?;
        let chunk_len = chunk.len();
        let relative_idx = idx - total_len;

        debug_assert_eq!(chunk_len, 0); // The length of a chunk should never be 0

        // match relative_idx

        if chunk_len == 1 {
            // Here we just remove the element if it is the only element in the array.
            self.refs.remove(chunk_idx);
            return Some(());
        }

        todo!();
        */

    }

    pub fn remove_range(&mut self, range: impl RangeBounds<usize>) -> Option<()> {

        let range_start = match range.start_bound() {
            Bound::Included(n) => *n,
            Bound::Excluded(n) => *n + 1,
            Bound::Unbounded => 0,
        };

        let range_end = match range.end_bound() {
            Bound::Included(n) => *n,
            Bound::Excluded(n) => *n - 1,
            Bound::Unbounded => 0,
        };

        let starting_chunk = self.get_chunk(chunk_idx)?;

        // let chunk = self.get_chunk(chunk_idx)?;
        let starting_chunk_len = starting_chunk.len();
        let relative_idx = idx - total_len;

        debug_assert_eq!(chunk_len, 0); // The length of a chunk should never be 0

        // match relative_idx

        if chunk_len == 1 {
            // Here we just remove the element if it is the only element in the array.
            self.refs.remove(chunk_idx);
            return Some(());
        }

        todo!();

    }

}

impl <T> From<Pin<Arc<[T]>>> for Sevec<T> {

    fn from(value: Pin<Arc<[T]>>) -> Self {
        // Gets the length
        let value_len = value.len();
        let ptr = ptr::slice_from_raw_parts(value.as_ptr(), value_len);
        return Self {
            data: vec![ value, ],
            refs: vec![ ptr,   ],
        };

    }
}

impl <T: std::fmt::Debug> std::fmt::Debug for Sevec<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        // Writes open bracket
        write!(f, "[")?;

        // Writes the inner data.
        for (i, ref_ptr) in self.refs.iter().enumerate() {

            // Goes through each slice
            let ref_slice = unsafe { &**ref_ptr };
            for (j, entry) in ref_slice.iter().enumerate() {
                // Writes value
                if i != 0 || j != 0 { write!(f, ", ")?; }
                write!(f, "{:?}", entry)?;
            }

        }

        // Writes closing bracket
        write!(f, "]")?;

        return Ok(());
    }
}

#[cfg(test)]
mod tests {


    use super::*;

    #[test]
    fn general_tests() {

        let mut v = Sevec::new();
        v.push("Hello There!");

        let vec = vec![
            "Hello There!",
        ];

        assert_eq!(format!("{:?}", vec), format!("{:?}", v));

    }

    #[test]
    fn test_getting_data() {

        let mut data = Sevec::new();

        // With Data Checks
        data.push(0);
        data.push(1);
        data.push(2);

        assert_eq!(data.get(0), Some(&0));
        assert_eq!(data.get(1), Some(&1));
        assert_eq!(data.get(2), Some(&2));
        assert_eq!(data.get(3), None);

        // Pushes new data
        data.push_arc_slice(Pin::new(vec![1, 2, 3].into()));

        // Previous Values
        assert_eq!(data.get(0), Some(&0));
        assert_eq!(data.get(1), Some(&1));
        assert_eq!(data.get(2), Some(&2));

        // Extended values
        assert_eq!(data.get(3), Some(&1));
        assert_eq!(data.get(4), Some(&2));
        assert_eq!(data.get(5), Some(&3));
        assert_eq!(data.get(6), None);

        // Data removed check
        data.remove(1).unwrap();

        assert_eq!(data.get(0), Some(&0));
        assert_eq!(data.get(1), Some(&2));
        assert_eq!(data.get(2), None);

        // Data removed check (from segmented section).
        todo!();
        // data.remove(1).unwrap();



    }

}
