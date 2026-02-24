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

    /// Gets the length of the inner data.
    /// This function is actually O(n) because we don't store the length as part of our structure.
    /// This may be done externally to improve performance however that is up to the implementor.
    pub fn len(&self) -> usize {
        return self.refs.iter()
            .map(|v| v.len())
            .sum()
            ;
    }

    // Inserts a new slice at a given chunk position.
    pub fn insert_arc_slice_to_chunk_pos(&mut self, chunk_index: usize, value: Pin<Arc<[T]>>) -> () {
        // Gets the reference
        let data_inner_ref = ptr::slice_from_raw_parts(value.as_ptr(), value.len());
        // Adds the data.
        self.data.push(value);
        // Adds the reference.
        self.refs.insert(chunk_index, data_inner_ref);
        return ();
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

    /*
    /// Calculates the estimated size of the sevec.
    fn size_estimation(&self) -> usize {

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
    */

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
    pub fn get_chunk_and_length_from_idx(&self, idx: usize) -> Option<(usize, usize)> {

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

        let (chunk_idx, total_len) = self.get_chunk_and_length_from_idx(idx)?;
        let chunk = self.get_chunk(chunk_idx)?;

        // Gets the result
        let final_idx = idx - total_len;
        let res = chunk.get(final_idx)?;
        return Some(res);

    }

    /// Removes the element at a given index.
    pub fn remove(&mut self, idx: usize) -> Option<()> {
        return self.remove_range(idx..=idx);
    }

    /// Removes all elements within the specified range.
    /// Without explicit bounds, this function will either start at the start or go until the end
    /// of all the data.
    ///
    /// ```rust
    /// # use sevec::Sevec;
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3, 4].into();
    /// sevec.remove_range(1..=2);
    /// assert_eq!(&format!("{:?}", sevec), "[1, 4]");
    /// ```
    pub fn remove_range(&mut self, range: impl RangeBounds<usize>) -> Option<()> {

        let range_start = match range.start_bound() {
            Bound::Included(&n) => n,
            Bound::Excluded(&n) => n + 1,
            Bound::Unbounded => 0,
        };

        let (starting_chunk_idx, starting_cumu_len) = self.get_chunk_and_length_from_idx(range_start)?;
        // This is the index of the start of the bounds within the start chunk
        let starting_chunk_rel_idx = range_start - starting_cumu_len;

        let range_end = match range.end_bound() {
            Bound::Unbounded => {
                // If the relative index is the start of a chunk.
                if starting_chunk_rel_idx == 0 {
                    // We remove the starting chunk and everything after it.
                    unsafe { self.refs.set_len(starting_chunk_idx); };
                    return Some(());
                }

                let start_mut = self.refs.get_mut(starting_chunk_idx)?;
                // Gets new location.
                *start_mut = ptr::slice_from_raw_parts_mut(start_mut.addr() as *mut _, starting_chunk_rel_idx);
                // Updates length
                unsafe { self.refs.set_len(starting_chunk_idx + 1); };
                return Some(());
            }

            Bound::Included(&n) => n,
            Bound::Excluded(&n) => n.checked_sub(1)?, // if n == 0
        };

        // This could be implemented a bit better considering we know starting_chunk_idx and
        // starting_cumu_len
        // Slight performance like this isn't a concern right now but should be considered in the
        // future. (TODO!)
        // We should have a function that works like
        // "get_chunk_and_length_from_idx_and_other_idx_and_len".
        // Very long name but this is okay because it would mainly be used internally (though
        // exposed externally).
        let (ending_chunk_idx, ending_cumu_len) = self.get_chunk_and_length_from_idx(range_end)?;
        let ending_chunk_rel_idx = range_end - ending_cumu_len;

        // Unsure if this is needed considering the behavior of range as well as how this will be
        // handled in the removal code.
        // // Bounds checking and such.
        // if ending_chunk_idx > starting_chunk_idx { return None; }

        // This unwrap shouldn't really ever fail.
        // This is between two indexes which are known good (or supposedly are).
        // If this fails then there are some serious problems with the state of the code.
        // Gets the updated first chunk.
        let mut starting_chunk = *self.refs.get(starting_chunk_idx).unwrap();
        starting_chunk = ptr::slice_from_raw_parts(starting_chunk.addr() as *const _, starting_chunk_rel_idx);

        // Gets the updated end chunk
        let mut ending_chunk = *self.refs.get(ending_chunk_idx).unwrap();
        ending_chunk = ptr::slice_from_raw_parts(
            (ending_chunk.addr() + (ending_chunk_rel_idx + 1) * size_of::<T>()) as *const _,
            ending_chunk.len() - (ending_chunk_rel_idx + 1)
        );

        // // Creates new array.
        // let mut out_refs = Vec::with_capacity(self.refs.len());

        let mut running_length = starting_chunk_idx;

        // Reserves enough room for one more element (we may add one more via splitting).
        // We do this just to make the `len` one more and to re-allocate for extra capacity.
        self.refs.push(ptr::slice_from_raw_parts(0 as *const _, 0));

        // Adds the chunks
        if starting_chunk.len() != 0 {
            self.refs[running_length] = starting_chunk;
            running_length += 1;
        }
        if ending_chunk.len() != 0 {
            self.refs[running_length] = ending_chunk;
            running_length += 1;
        }

        // This might be able to be replaced with a [`ptr::copy`] call however, in many cases this
        // might just be shifting one element at a time where the speedups may be very little.
        for i in (ending_chunk_idx + 1)..self.refs.len() {
            self.refs[running_length] = self.refs[i];
            running_length += 1;
        }

        // Update the length
        unsafe { self.refs.set_len(running_length); };

        return Some(());

    }

}

impl <T: Unpin> Sevec<T> {

    /// Adds a value to the end of the array.
    pub fn push(&mut self, value: T) -> () {

        // Creates new ptr.
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
        return;

    }

}

impl <T: Unpin + Clone + Sized> Sevec<T> {

    /// Copies and inserts a given slice.
    pub fn push_slice(&mut self, value: &[T]) -> () {

        let arc_ptr = Arc::<[T]>::new_uninit_slice(value.len());
        let mut arc_ptr = unsafe { arc_ptr.assume_init() };

        let arc_ptr_mut = Arc::get_mut(&mut arc_ptr).unwrap();
        let arc_ptr_mut_raw = arc_ptr_mut.as_mut_ptr();

        // Copies the data.
        unsafe {
            ptr::copy_nonoverlapping(value.as_ptr(), arc_ptr_mut_raw, value.len());
        }

        let arc_ptr = Pin::new(arc_ptr);

        self.push_arc_slice(arc_ptr);

        return;

    }

}

impl <T: std::fmt::Debug> std::fmt::Debug for Sevec<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        // Writes open bracket
        write!(f, "[")?;

        // Flag for if the item being written is the first.
        let mut first = true;

        // Writes the inner data.
        for ref_ptr in self.refs.iter() {

            // Goes through each slice
            let ref_slice = unsafe { &**ref_ptr };
            for entry in ref_slice.iter() {
                // Checks if value is first item
                match first {
                    true  => { first = false; },
                    false => write!(f, ", ")?,
                }
                // Writes value
                write!(f, "{:?}", entry)?;
            }

        }

        // Writes closing bracket
        write!(f, "]")?;

        return Ok(());
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

impl <T: Unpin + Clone + Sized> From<&[T]> for Sevec<T> {
    fn from(value: &[T]) -> Self {
        let mut data = Self::new();
        data.push_slice(value);
        return data;
    }
}

impl <T: Unpin + Clone + Sized> From<Vec<T>> for Sevec<T> {
    fn from(value: Vec<T>) -> Self {
        let slice = value.as_slice();
        return slice.into();
    }
}

impl <T: Clone> Into<Vec<T>> for &Sevec<T> {
    fn into(self) -> Vec<T> {

        // Technically self.len is O(n) but the O(n) is pretty fast.
        // I imagine this is likely worth not re-allocating over and over but that is just an
        // assumption.
        let out_len = self.len();
        let mut new_vec = Vec::<T>::with_capacity(out_len);
        let new_ptr_addr = new_vec.as_mut_ptr();
        let mut length_sum = 0;

        for chunk in &self.refs {
            // Copies Data
            unsafe {
                ptr::copy_nonoverlapping(
                    chunk.addr() as *const T,
                    new_ptr_addr.add(length_sum),
                    chunk.len()
                );
            };
            // Updates last byte
            length_sum += chunk.len();
        }

        // Updates vec length
        unsafe{ new_vec.set_len(out_len); }

        return new_vec;

    }
}

impl <T: Clone> Into<Vec<T>> for Sevec<T> {
    fn into(self) -> Vec<T> {
        return (&self).into();
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_display() {
        let mut v = Sevec::new();
        v.push("Hello There!");
        let vec = vec![
            "Hello There!",
        ];
        assert_eq!(format!("{:?}", vec), format!("{:?}", v));
        v.insert_arc_slice_to_chunk_pos(0, Pin::new(vec!["Hello", "H"].into()));
        v.remove_range(0..2).unwrap();
        assert_eq!(format!("{:?}", vec), format!("{:?}", v));
    }

    #[test]
    fn test_remove_basic() {
        let mut data = Sevec::new();
        data.push(1);
        data.push(2);
        data.push(3);
        data.push(4);
        data.remove_range(1..=2).unwrap();
        assert_eq!(data.len(), 2);
        assert_eq!(data.get(0), Some(&1));
        assert_eq!(data.get(1), Some(&4));
    }

    #[test]
    fn test_remove_out_of_range() {
        let mut data = Sevec::new();
        data.push(1);
        data.push(2);
        data.push(3);
        data.push(4);

        assert!(data.remove_range(0..5).is_none());
        assert!(data.remove_range(0..100).is_none());

        let out_data: Vec<_> = data.clone().into();
        assert_eq!(out_data, vec![1, 2, 3, 4]);

        let res = data.remove_range(0..4);
        assert!(res.is_some());

        assert_eq!(data.len(), 0);
    }

    #[test]
    fn test_remove_everything() {
        let mut data = Sevec::new();
        data.push(1);
        data.push(2);
        data.push(3);
        data.push(4);
        data.remove_range(0..data.len()).unwrap();
        assert_eq!(data.len(), 0);
        // This implies that no empty pointers exist.
        // Not really a required attribute but is nice to have.
        // assert_eq!(data.refs.len(), 0);
    }

    #[test]
    fn test_remove_slices() {

        let mut data = Sevec::new();
        data.push_arc_slice(Pin::new(vec![5, 2, 1].into()));
        data.push(1);
        data.push(2);
        data.push(3);
        data.push(4);
        data.push_arc_slice(Pin::new(vec![1, 2, 3].into()));

        assert_eq!(format!("{:?}", data), format!("{:?}",
            vec![5, 2, 1, 1, 2, 3, 4, 1, 2, 3]
        ));

        data.remove_range(1..9).unwrap(); // Remove across two slices
        assert_eq!(format!("{:?}", data), format!("{:?}", vec![5, 3]));
        data.remove_range(1..).unwrap(); // Unbounded remove
        assert_eq!(format!("{:?}", data), format!("{:?}", vec![5]));
        data.push(3);
        data.remove_range(..data.len()).unwrap(); // Unbounded remove from front
        assert_eq!(format!("{:?}", data), format!("{:?}", Vec::<i32>::new()));
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
        assert_eq!(data.get(5), None); // After removal the end would have moved.

        data.remove_range(..).unwrap();
        assert_eq!(data.len(), 0);

    }

    #[test]
    fn test_push_slice() {

        let mut data: Sevec<u32> = vec![1, 23, 3, 3].into();
        data.push_slice(&[1, 2, 3, 4]);

        let reference_vec: Vec<_> = data.into();
        assert_eq!(
            reference_vec,
            vec![1, 23, 3, 3, 1, 2, 3, 4]
        );

    }

    #[test]
    fn test_remove_out_of_range() {
        todo!();
    }

}
