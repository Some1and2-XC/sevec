#![doc = include_str!("../README.md")]

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
    /// This may be done externally to improve performance, however, that is up to library consumers.
    /// ```rust
    /// # use sevec::Sevec;
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3].into();
    /// assert_eq!(sevec.len(), 3);
    /// ```
    pub fn len(&self) -> usize {
        return self.refs.iter()
            .map(|v| v.len())
            .sum()
            ;
    }

    /// Gets a reference to some data.
    /// ```rust
    /// # use sevec::Sevec;
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3].into();
    /// assert_eq!(sevec.get(0), Some(&1));
    /// assert_eq!(sevec.get(1), Some(&2));
    /// assert_eq!(sevec.get(2), Some(&3));
    /// assert_eq!(sevec.get(3), None);
    /// ```
    pub fn get(&self, idx: usize) -> Option<&T> {
        let (chunk_idx, total_len) = self.get_chunk_and_length_from_idx(idx)?;
        let chunk = self.get_chunk(chunk_idx)?;
        // Gets the result
        let final_idx = idx - total_len;
        let res = chunk.get(final_idx)?;
        return Some(res);
    }

}

impl <T: Unpin> Sevec<T> {

    /// Adds a value to the end of the array.
    /// ```rust
    /// # use sevec::Sevec;
    /// let mut sevec = Sevec::new();
    /// sevec.push(1);
    /// assert_eq!(sevec.to_string(), "[1]");
    /// assert_eq!(sevec.len(), 1);
    /// ```
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
    /// ```rust
    /// # use sevec::Sevec;
    /// let mut sevec = Sevec::new();
    /// sevec.push_slice(&[1, 2, 3, 4]);
    /// assert_eq!(sevec.to_string(), "[1, 2, 3, 4]");
    /// assert_eq!(sevec.len(), 4);
    /// ```
    pub fn push_slice(&mut self, value: &[T]) -> () {

        let arc_ptr = Arc::<[T]>::new_uninit_slice(value.len());
        let mut arc_ptr = unsafe { arc_ptr.assume_init() };
        let arc_ptr_mut = Arc::get_mut(&mut arc_ptr).unwrap();

        // Copies the data.
        unsafe {
            ptr::copy_nonoverlapping(value.as_ptr(), arc_ptr_mut.as_mut_ptr(), value.len());
        }

        self.push_arc_slice(Pin::new(arc_ptr));

        return;

    }

    /// Creates a copy of the slice data and inserts it at a specified location.
    /// This method uses [`Self::insert_arc_slice`] and if able, this method should be prefered
    /// with direct data writing to the underlying [`Arc`] data structure.
    /// ```rust
    /// # use sevec::Sevec;
    /// // Initializes the data
    /// let mut sevec = vec![1, 2, 3, 4].into();
    ///
    /// // Inserts numbers between the 1 and 2
    /// sevec.insert_slice(1, &[100, 200]).unwrap();
    ///
    /// // Displays the result
    /// assert_eq!(sevec.to_string(), "[1, 100, 200, 2, 3, 4]");
    /// ```
    pub fn insert_slice(&mut self, idx: usize, value: &[T]) -> Option<()> {

        let arc_ptr = Arc::<[T]>::new_uninit_slice(value.len());
        let mut arc_ptr = unsafe { arc_ptr.assume_init() };
        let arc_ptr_mut = Arc::get_mut(&mut arc_ptr).unwrap();

        // Copies the data.
        unsafe {
            ptr::copy_nonoverlapping(value.as_ptr(), arc_ptr_mut.as_mut_ptr(), value.len());
        }

        return self.insert_arc_slice(idx, Pin::new(arc_ptr));
    }

}

impl <T: Clone> Into<Vec<T>> for Sevec<T> {
    fn into(self) -> Vec<T> {
        return (&self).into();
    }
}

impl <T> Sevec<T> {

    /// Adds a new slice.
    /// This method is particularly well suited to situations where direct writing to the inner
    /// value of an [`Arc`] pointer is available. This method moves the [`Arc`] pointer without
    /// copying.
    ///
    /// The following example is relatively long on account of it showcasing how to create an
    /// uninitialized [`Arc`] pointer with the purpose of minimizing copies.
    /// ```rust
    /// # use sevec::Sevec;
    /// # use std::{mem::MaybeUninit, pin::Pin, sync::Arc};
    /// let mut sevec: Sevec<u32> = Sevec::new();
    ///
    /// let data_len = 6; // Example array size
    ///
    /// // Creates the pointer uninitialized to avoid zeroing.
    /// let mut data = {
    ///     let ptr = Arc::<[u32]>::new_uninit_slice(data_len);
    ///     unsafe { ptr.assume_init() }
    /// };
    ///
    /// // Here we get mutable access to the data. We call unwrap without
    /// // worry because we know we have exclusive access to the pointer.
    /// let data_mut = Arc::get_mut(&mut data).unwrap();
    ///
    /// // Writing the data directly to the [`Arc`] ptr.
    /// for i in 0..data_mut.len() {
    ///     data_mut[i] = i as u32;
    /// }
    ///
    /// // Calling this function to push the data into the [`Sevec`].
    /// sevec.push_arc_slice(Pin::new(data));
    ///
    /// // We can see the values we wrote directly into the [`Arc`] get displayed.
    /// assert_eq!(sevec.to_string(), "[0, 1, 2, 3, 4, 5]");
    /// ```
    pub fn push_arc_slice(&mut self, value: Pin<Arc<[T]>>) -> () {
        // Gets the reference
        let data_inner_ref = ptr::slice_from_raw_parts(value.as_ptr(), value.len());
        // Adds the data.
        self.data.push(value);
        // Adds the reference.
        self.refs.push(data_inner_ref);
        return ();
    }

    /// Inserts a slice into a specified location in the [`Sevec`].
    /// ```rust
    /// # use sevec::Sevec;
    /// # use std::{pin::Pin, sync::Arc};
    /// // Creates array
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3].into();
    /// // Creates data
    /// let data: Arc<[u32]> = vec![4, 5, 6].into_boxed_slice().into();
    ///
    /// // Inserts data between the `1` and `2`.
    /// sevec.insert_arc_slice(1, Pin::new(data));
    ///
    /// // Shows result
    /// assert_eq!(&sevec.to_string(), "[1, 4, 5, 6, 2, 3]");
    /// ```
    pub fn insert_arc_slice(&mut self, idx: usize, value: Pin<Arc<[T]>>) -> Option<()> {
        // Gets slice
        let slice = ptr::slice_from_raw_parts(value.as_ptr(), value.len());
        // Tries to insert.
        unsafe { self.insert_raw_slice(idx, slice) }?;
        // Inserts Data only if the slice was added.
        self.data.push(value);
        // Returns result.
        return Some(());
    }

    /// Adds a slice to the array.
    /// This should only be done with a slice which has a lifetime associated with the lifetime of
    /// this object, in particular, either a slice refering to something with a static lifetime or
    /// containing data within the data of this array is intended.
    ///
    /// This function is intended to be used in cases where repeated data is going to be added to
    /// the list and therefore adding the data repeatedly to the inner data stores isn't needed.
    /// ```rust
    /// # use sevec::Sevec;
    /// # use std::{ptr, pin::Pin, sync::Arc};
    /// // Creates array
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3].into();
    /// // Creates data
    /// static DATA: &'static [u32; 3] = &[4, 5, 6];
    ///
    /// // Adds the data between the `1` and `2`
    /// unsafe {
    ///     sevec.insert_raw_slice(1, ptr::slice_from_raw_parts(DATA.as_ptr(), DATA.len()));
    /// };
    ///
    /// // Checks result
    /// assert_eq!(&sevec.to_string(), "[1, 4, 5, 6, 2, 3]");
    /// ```
    pub unsafe fn insert_raw_slice(&mut self, idx: usize, slice: *const [T]) -> Option<()> {

        let (mut write_idx, left, right) = match self.get_chunk_and_length_from_idx(idx) {
            Some((chunk_idx, prev_sum))=> {

                let chunk = *self.refs.get(chunk_idx)?;
                let offset = idx - prev_sum;

                // Creates left and right sides
                let left = ptr::slice_from_raw_parts(chunk as *const T, offset);
                let right = ptr::slice_from_raw_parts((chunk as *const T).add(offset), chunk.len() - offset);

                (chunk_idx, left, right)
            },
            None => {
                // If we didn't find the chunk but the chunk index was 0, we continue.
                if idx != 0 {
                    return None;
                }
                (0, ptr::slice_from_raw_parts(ptr::null(), 0), ptr::slice_from_raw_parts(ptr::null(), 0))
            },
        };

        // Inserts values

        // We write this to the left of the chunk (if there was a chunk)
        if left.len() != 0 {
            self.refs.insert(write_idx, left);
            write_idx += 1;
        }

        // We write the actual data overtop of the original chunk each time.
        if slice.len() != 0 {
            if let Some(data) = self.refs.get_mut(write_idx) {
                *data = slice;
            }
            else {
                self.refs.insert(write_idx, slice);
            }

            // self.refs.insert(write_idx, slice);
            write_idx += 1;
        }

        if right.len() != 0 {
            self.refs.insert(write_idx, right);
            // write_idx += 1;
        }

        return Some(());

    }

    /// Removes the element at a given index.
    /// ```rust
    /// # use sevec::Sevec;
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3].into();
    /// sevec.remove(1); // Removes the middle `2`
    /// assert_eq!(sevec.to_string(), "[1, 3]");
    /// ```
    pub fn remove(&mut self, idx: usize) -> Option<()> {
        return self.remove_range(idx..=idx);
    }

    /// Removes all elements within the specified range.
    /// Without explicit bounds, this function will either start at the start or go until the end
    /// of all the data.
    ///
    /// This function is a more ergonomic way of calling [`Self::remove_between_start_and_end`],
    /// that funciton is also available if the [`RangeBounds`] handling overhead is unwanted.
    ///
    /// ```rust
    /// # use sevec::Sevec;
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3, 4].into();
    /// sevec.remove_range(1..=2); // Removes both `2` and `3`
    /// assert_eq!(sevec.to_string(), "[1, 4]");
    /// assert_eq!(sevec.len(), 2);
    /// ```
    pub fn remove_range(&mut self, range: impl RangeBounds<usize>) -> Option<()> {

        let range_start = match range.start_bound() {
            Bound::Included(&n) => n,
            Bound::Excluded(&n) => n + 1,
            Bound::Unbounded => 0,
        };

        let range_end = match range.end_bound() {
            Bound::Unbounded => {

                let (starting_chunk_idx, starting_cumu_len) = self.get_chunk_and_length_from_idx(range_start)?;
                // This is the index of the start of the bounds within the start chunk
                let starting_chunk_rel_idx = range_start - starting_cumu_len;

                // If the relative index is the start of a chunk.
                if starting_chunk_rel_idx == 0 {
                    // We remove the starting chunk and everything after it.
                    unsafe { self.refs.set_len(starting_chunk_idx); };
                    return Some(());
                }

                let start_mut = self.refs.get_mut(starting_chunk_idx)?;
                // Gets new location.
                *start_mut = ptr::slice_from_raw_parts_mut(*start_mut as *mut _, starting_chunk_rel_idx);
                // Updates length
                unsafe { self.refs.set_len(starting_chunk_idx + 1); };
                return Some(());
            }

            Bound::Included(&n) => n,
            Bound::Excluded(&n) => n.checked_sub(1)?, // if n == 0
        };

        return self.remove_between_start_and_end(range_start, range_end);

    }

    /// Removes all elements within the specified range.
    /// Both the start and end values are inclusive.
    ///
    /// For an ergonomic wrapper around this method, use [`Self::remove_range`].
    ///
    /// ```rust
    /// # use sevec::Sevec;
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3, 4].into();
    /// // Removes everything between index 1 and 2 (numbers `2` and `3`)
    /// sevec.remove_between_start_and_end(1, 2);
    /// assert_eq!(sevec.to_string(), "[1, 4]");
    /// assert_eq!(sevec.len(), 2);
    /// ```
    pub fn remove_between_start_and_end(&mut self, range_start: usize, range_end: usize) -> Option<()> {

        let (starting_chunk_idx, starting_cumu_len) = self.get_chunk_and_length_from_idx(range_start)?;
        // This is the index of the start of the bounds within the start chunk
        let starting_chunk_rel_idx = range_start - starting_cumu_len;

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

        // This unwrap shouldn't really ever fail.
        // This is between two indexes which are known good (or supposedly are).
        // If this fails then there are some serious problems with the state of the code.
        // Gets the updated first chunk.
        let mut starting_chunk = *self.refs.get(starting_chunk_idx).unwrap();
        starting_chunk = ptr::slice_from_raw_parts(starting_chunk as *const _, starting_chunk_rel_idx);

        // Gets the updated end chunk
        let mut ending_chunk = *self.refs.get(ending_chunk_idx).unwrap();
        ending_chunk = ptr::slice_from_raw_parts(
            unsafe {
                (ending_chunk as *const T).add(ending_chunk_rel_idx + 1)
            },
            ending_chunk.len() - (ending_chunk_rel_idx + 1)
        );

        // // Creates new array.
        // let mut out_refs = Vec::with_capacity(self.refs.len());

        let mut running_length = starting_chunk_idx;

        // Reserves enough room for one more element (we may add one more via splitting).
        // We do this just to make the `len` one more and to re-allocate for extra capacity.
        self.refs.push(ptr::slice_from_raw_parts(ptr::null(), 0));

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
        // We subtract 1 from the length because of the padded value added earlier.
        for i in (ending_chunk_idx + 1)..(self.refs.len() - 1) {
            self.refs[running_length] = self.refs[i];
            running_length += 1;
        }

        // Update the length
        unsafe { self.refs.set_len(running_length); };

        return Some(());

    }

    /// Gets a specified chunk as a slice.
    /// Note, this is the underlying chunk, not the actual data at a given index.
    /// ```rust
    /// # use sevec::Sevec;
    /// // Initializes the array
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3].into();
    ///
    /// // Pushes more data
    /// sevec.push_slice(&[4, 5, 6]);
    ///
    /// // Gets the initial data
    /// let data = sevec.get_chunk(0).unwrap();
    ///
    /// // Checks result
    /// assert_eq!(&data, &[1, 2, 3]); // We got the first allocation only.
    /// ```
    pub fn get_chunk(&self, chunk: usize) -> Option<&[T]> {
        let chunk_ptr = self.refs.get(chunk)?;
        return unsafe { chunk_ptr.as_ref() };
    }

    /// Gets a specified value from both a chunk index and a chunk sub index.
    /// ```rust
    /// # use sevec::Sevec;
    /// // Initializes the array
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3].into();
    /// assert_eq!(*sevec.get_from_chunk_and_idx(0, 0).unwrap(), 1); // Gets the first item
    /// assert_eq!(sevec.get_from_chunk_and_idx(1, 0), None); // Doesn't Exist yet
    ///
    /// // Pushes more data
    /// sevec.push_slice(&[4, 5, 6]);
    /// // Gets the first item from the second allocation
    /// assert_eq!(*sevec.get_from_chunk_and_idx(1, 0).unwrap(), 4);
    /// ```
    pub fn get_from_chunk_and_idx(&self, chunk: usize, idx: usize) -> Option<&T> {
        let chunk_slice = self.get_chunk(chunk)?;
        return chunk_slice.get(idx);
    }

    /// Gets the chunk index of a specified input index.
    /// The first value in the result is the chunk index.
    /// The second value is the sum of all the previous lengths up until the specified chunk.
    /// ```rust
    /// # use sevec::Sevec;
    /// // Initializes the array
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3].into();
    ///
    /// // Adds more data
    /// sevec.push_slice(&[4, 5, 6]);
    /// assert_eq!(&sevec.to_string(), "[1, 2, 3, 4, 5, 6]");
    ///
    /// // Gets the location of index 3
    /// let (chunk_idx, prev_sum) = sevec.get_chunk_and_length_from_idx(3).unwrap();
    ///
    /// // The chunk index is 1
    /// assert_eq!(chunk_idx, 1);
    ///
    /// // prev_sum is the amount of data that appeared before the chunk that was discovered.
    /// assert_eq!(prev_sum, 3); // Because the first allocation had 3 items.
    /// ```
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

    /// Inserts a new slice at a given chunk position.
    /// This method is especially useful if data can be written directly into the [`Arc<[T]>`].
    /// An example of how this can be done can be found in the docs of [`Self::push_arc_slice`].
    /// It is important to note that this method works on the internal chunk position, rather than
    /// the index of the array.
    ///
    /// It is also likely that the method wanted in this case is actually
    /// [`Self::insert_arc_slice`].
    /// ```rust
    /// # use sevec::Sevec;
    /// # use std::{pin::Pin, sync::Arc};
    /// // Creating a sevec
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3].into();
    ///
    /// // Creating new data
    /// let slice_data = Arc::new([4, 5, 6]);
    ///
    /// // Adding the data before anything else.
    /// sevec.insert_arc_slice_to_chunk_pos(0, Pin::new(slice_data));
    ///
    /// // Checking the result.
    /// assert_eq!(&sevec.to_string(), "[4, 5, 6, 1, 2, 3]");
    /// ```
    pub fn insert_arc_slice_to_chunk_pos(&mut self, chunk_index: usize, value: Pin<Arc<[T]>>) -> () {
        // Gets the reference
        let data_inner_ref = ptr::slice_from_raw_parts(value.as_ptr(), value.len());
        // Adds the data.
        self.data.push(value);
        // Adds the reference.
        self.refs.insert(chunk_index, data_inner_ref);
        return ();
    }

}

impl <T: std::fmt::Debug> std::fmt::Debug for Sevec<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        // Writes open bracket
        write!(f, "[")?;

        // Flag for if the item being written is the first.
        let mut first = true;

        // Writes the inner data.
        for &ref_ptr in self.refs.iter() {

            // Goes through each slice
            let ref_slice = unsafe { &*ref_ptr };

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

impl <T: std::fmt::Debug> std::string::ToString for Sevec<T> {
    fn to_string(&self) -> String {
        return format!("{:?}", self);
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
                    *chunk as *const T,
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
    fn test_insert_raw_slice() {
        let mut sevec = Sevec::new();
        sevec.push_slice(&[1, 2, 3]);
        let data = &[4, 5, 6];
        let data_ptr = ptr::slice_from_raw_parts(data.as_ptr(), data.len());
        unsafe {
            sevec.insert_raw_slice(1, data_ptr)
        }.unwrap();
        assert_eq!(sevec.to_string(), "[1, 4, 5, 6, 2, 3]");

        sevec.remove_range(..);

        assert!(unsafe { sevec.insert_raw_slice(1, &[1, 2, 3]) }.is_none());
        assert_eq!(sevec.len(), 0);
        assert_eq!(sevec.refs.len(), 0);

        unsafe { sevec.insert_raw_slice(0, data_ptr) }.unwrap();
        assert_eq!(&sevec.to_string(), "[4, 5, 6]");

        unsafe { sevec.insert_raw_slice(0, data_ptr) }.unwrap();
        assert_eq!(&sevec.to_string(), "[4, 5, 6, 4, 5, 6]");

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
    fn test_remove_other_range() {

        let mut data: Sevec<_> = vec![1, 2, 3, 4].into();
        let res = data.remove_range(1..=2);
        assert!(res.is_some());

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
        assert_eq!(data.refs.len(), 0);
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

}
