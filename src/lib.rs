#![doc = include_str!("../README.md")]

use std::{cmp::Ordering, ops::{Bound, RangeBounds}, ptr, sync::Arc};

#[derive(Clone)]
pub struct Sevec<T> {
    /// The underlying data.
    data: Vec<Arc<[T]>>,
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

        // Pushes the value
        // In theory, we could just hard-code adding the slice length of 1 but I don't really
        // think this would make much of a difference.
        self.push_arc_slice(arc_ptr);
        return;

    }

    /// Joins two [`Sevec<T>`]'s together.
    /// The `other` value gets added onto the end.
    /// ```rust
    /// # use sevec::Sevec;
    /// let mut sevec: Sevec<u8> = vec![1, 2, 3].into();
    /// sevec.extend(vec![4, 5, 6].into());
    /// assert_eq!(sevec.to_string(), "[1, 2, 3, 4, 5, 6]");
    /// ```
    pub fn extend(&mut self, other: Self) -> () {
        let Self { data, refs} = other;
        // Pushes data
        for v in data.into_iter() {
            self.data.push(v);
        }
        // Pushes refs
        for v in refs.into_iter() {
            self.refs.push(v);
        }
    }

    /// Gets ordinary references to the underlying data.
    /// This is done by memory transmutation and may be compiled out.
    #[inline]
    pub fn get_refs(&self) -> &[&[T]] {

        // SAFETY: the memory layout of a const *[T] should be the same as a &[T].
        // Therefore, the operation created the proper datatype.
        // For the lifetime of the returned data, the underlying data should last
        // no longer than the underlying data. Modification of the underlying data should
        // also not be possible while the created reference exists.
        return unsafe {
            // TODO: I'm not sure that the existing "SAFTEY" comment is good enough at this point.
            // Using [`std::mem::transmute`] is very unsafe and variance should be explored
            // further. That being said --- it's probably fine.
            std::mem::transmute(&self.refs as &[_])
        };

    }

    /// This method inserts a slice into the specified location so long as the slice has a 'static
    /// lifetime. Similar to the [`Sevec::insert_slice`] method with the primary difference is that
    /// the slice doesn't get copied into an internal [`Arc`] allocation.
    ///
    /// Because the slice data doesn't need to be copied (because it has a 'static lifetime), T
    /// doesn't need to implement [`Clone`] either.
    ///
    /// ```rust
    /// # use sevec::Sevec;
    /// // Creates the `Sevec<T>`
    /// let mut sevec: Sevec<_> = vec![1, 6, 7, 8].into();
    /// assert_eq!(sevec.to_string(), "[1, 6, 7, 8]");
    ///
    /// // Inserts data into index 1
    /// const DATA: &[usize] = &[2, 3, 4, 5];
    /// sevec.insert_static_slice(1, DATA).unwrap();
    ///
    /// // Validates the result
    /// assert_eq!(sevec.to_string(), "[1, 2, 3, 4, 5, 6, 7, 8]");
    /// assert_eq!(sevec.len(), 8);
    /// ```
    pub fn insert_static_slice(&mut self, idx: usize, value: &'static [T]) -> Option<()> {
        // Creates raw slice from data.
        let slice = ptr::slice_from_raw_parts(value.as_ptr(), value.len());
        // Pushes slice into place.
        return unsafe {
            self.insert_raw_slice(idx, slice)
        };
    }

}

impl <T: Clone + Sized> Sevec<T> {

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

        self.push_arc_slice(arc_ptr);

        return;

    }

    /// Creates a copy of the slice data and inserts it at a specified location.
    /// This method uses [`Self::insert_arc_slice`] and if able, this method should be preferred
    /// with direct data writing to the underlying [`Arc`] data structure.
    /// ```rust
    /// # use sevec::Sevec;
    /// // Initializes the data
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3, 4].into();
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

        arc_ptr_mut.clone_from_slice(value);

        return self.insert_arc_slice(idx, arc_ptr);
    }

    /// Removes the end of a slice and copies to out buffer.
    /// ```rust
    /// # use sevec::Sevec;
    /// // Initializes the array
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3].into();
    ///
    /// // Slices data
    /// let res = sevec.remove_and_copy_slice_from_end(2).unwrap();
    ///
    /// // We are left with the first allocation only
    /// assert_eq!(&sevec.to_string(), "[1]");
    /// // The result is the last `amnt` elements (in this case 2 elements)
    /// assert_eq!(&*res, &[2, 3]);
    /// ```
    pub fn remove_and_copy_slice_from_end(&mut self, amnt: usize) -> Option<Arc<[T]>> {

        let mut amnt_sliced = 0;
        let mut current_idx = self.refs.len();

        let out_data = Arc::new_uninit_slice(amnt);
        let mut out_data: Arc<[T]> = unsafe { out_data.assume_init() };
        let out_data_mut = Arc::get_mut(&mut out_data).unwrap();

        let mut cur_ref_iter;

        loop {

            if current_idx == 0 {
                return None;
            }

            current_idx -= 1;

            let cur_ref = self.refs[current_idx];
            amnt_sliced += cur_ref.len();
            cur_ref_iter = unsafe { cur_ref.as_ref() }.unwrap();

            match amnt_sliced.cmp(&amnt) {
                Ordering::Less => {
                    for (i, v) in cur_ref_iter.iter().enumerate() {
                        out_data_mut[amnt - amnt_sliced + i] = v.clone();
                    }
                },
                Ordering::Equal => {
                    for (i, v) in cur_ref_iter.iter().enumerate() {
                        out_data_mut[amnt - amnt_sliced + i] = v.clone();
                    }
                    break;
                },
                Ordering::Greater => {
                    let diff = amnt_sliced - amnt;
                    for (i, v) in cur_ref_iter[diff..].iter().enumerate() {
                        out_data_mut[i] = v.clone();
                    }
                    break;
                },
            }

        }

        // Updates length
        // SAFETY: This is done in this way because we know we are just shrinking the vec.
        // Because `current_idx` is always less than `self.refs.len()`, out of range memory
        // is never referred to.
        let _ = unsafe { self.refs.set_len(current_idx) };

        // IF we cut one in half, we add back the start.
        if amnt_sliced > amnt {
            // If we need to cut one in half.
            let diff = amnt_sliced - amnt;
            self.refs.push(
                ptr::slice_from_raw_parts(cur_ref_iter.as_ptr(), diff)
            );
        }

        return Some(out_data);

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
    /// # use std::{mem::MaybeUninit, sync::Arc};
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
    /// sevec.push_arc_slice(data);
    ///
    /// // We can see the values we wrote directly into the [`Arc`] get displayed.
    /// assert_eq!(sevec.to_string(), "[0, 1, 2, 3, 4, 5]");
    /// ```
    pub fn push_arc_slice(&mut self, value: Arc<[T]>) -> () {
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
    /// # use std::sync::Arc;
    /// // Creates array
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3].into();
    /// // Creates data
    /// let data: Arc<[u32]> = vec![4, 5, 6].into_boxed_slice().into();
    ///
    /// // Inserts data between the `1` and `2`.
    /// sevec.insert_arc_slice(1, data);
    ///
    /// // Shows result
    /// assert_eq!(&sevec.to_string(), "[1, 4, 5, 6, 2, 3]");
    /// ```
    pub fn insert_arc_slice(&mut self, idx: usize, value: Arc<[T]>) -> Option<()> {
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
    /// this object, in particular, either a slice referring to something with a static lifetime or
    /// containing data within the data of this array is intended.
    ///
    /// If a slice has a static lifetime however, using [`Sevec::insert_static_slice`] is preferred
    /// over using this unsafe function.
    ///
    /// This function is intended to be used in cases where repeated data is going to be added to
    /// the list and therefore adding the data repeatedly to the inner data stores isn't needed.
    /// ```rust
    /// # use sevec::Sevec;
    /// # use std::{ptr, sync::Arc};
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

                // By this point, we have found a valid index so if we have a length of 0, we
                // skip everything and just say the operation would have been successful.
                if slice.len() == 0 {
                    return Some(());
                }

                let chunk = *self.refs.get(chunk_idx)?;
                let offset = idx - prev_sum;

                // Creates left and right sides
                let left = ptr::slice_from_raw_parts(chunk as *const T, offset);
                let right = ptr::slice_from_raw_parts((chunk as *const T).add(offset), chunk.len() - offset);

                (chunk_idx, left, right)
            },
            None => {

                // If we are at the last element, we insert anyway.
                // This is another O(n) operation, however, because this is the bad path anyway the
                // performance isn't too much of a concern.
                // Also, the compiler might do its thing and make this not exist in this way anyway.
                //
                // It is a shame that we can't get the total length out of
                // [`Sevec::get_chunk_and_length_from_idx`] but it really isn't worth the API
                // change to do convert the return type from `Option<(usize, usize)>` to something
                // like `(usize, Option<usize>)`.
                if idx == self.len() {
                    if slice.len() != 0 {
                        self.refs.push(slice);
                    }
                    return Some(());
                }
                else {
                    return None;
                }

            },
        };

        // Inserts values

        // We write this to the left of the chunk (if there was a chunk)
        if left.len() != 0 {
            self.refs.insert(write_idx, left);
            write_idx += 1;
        }

        // We write the actual data over-top of the original chunk each time.
        if let Some(data) = self.refs.get_mut(write_idx) {
            *data = slice;
        }
        else {
            self.refs.insert(write_idx, slice);
        }
        write_idx += 1;

        // Finally, we write the right side.
        if right.len() != 0 {
            self.refs.insert(write_idx, right);
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
    /// that function is also available if the [`RangeBounds`] handling overhead is unwanted.
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

        // If the range is backwards, return early.
        if range_start > range_end {
            return None;
        }

        // Gets chunks and chunk indexes

        let mut len_cumu = 0;
        let mut starting_chunk_idx = 0usize;
        // Initializes the chunks to null.
        // SAFETY: This is fine so long as neither actual chunk gets de-referenced, and we check
        // the length before adding each chunk to it's respective location.
        let mut starting_chunk = ptr::slice_from_raw_parts(ptr::null(), 0);
        let mut ending_chunk = ptr::slice_from_raw_parts(ptr::null(), 0);

        // Finds starting chunk
        for &ptr in self.get_refs() {
            let addr = ptr.as_ptr();
            let len = ptr.len();
            len_cumu += len;
            if len_cumu >= range_start {
                len_cumu -= len;
                starting_chunk = ptr::slice_from_raw_parts(addr, range_start - len_cumu);
                // if len_cumu - range_start >= 1 {
                // }

                break;
            }
            starting_chunk_idx += 1;
        }

        // Early return on not found.
        if starting_chunk_idx >= self.refs.len() {
            return None;
        }

        let mut ending_chunk_idx = starting_chunk_idx;
        // Finds the end chunk
        for &ptr in self.get_refs()[starting_chunk_idx..].iter() {
            let addr = ptr.as_ptr();
            let len = ptr.len();
            len_cumu += len;
            if len_cumu > range_end {
                len_cumu -= len;
                // let chunk_len = len_cumu - range_end;
                let chunk_len = range_end - len_cumu;
                ending_chunk = ptr::slice_from_raw_parts(
                    unsafe { addr.add(chunk_len + 1) },
                    len - (chunk_len + 1)
                );
                break;
            }
            ending_chunk_idx += 1;
        }

        // Early return on not found.
        if ending_chunk_idx >= self.refs.len() {
            return None;
        }

        // --- End of getting indexes and chunks ----

        // Case where we may be adding an entry.
        if starting_chunk_idx == ending_chunk_idx {

            match (starting_chunk.len(), ending_chunk.len()) {
                (0, 0) => {
                    self.refs.remove(starting_chunk_idx);
                },
                (0, _) => {
                    self.refs[starting_chunk_idx] = ending_chunk;
                },
                (_, 0) => {
                    self.refs[starting_chunk_idx] = starting_chunk;
                },
                (_, _) => {
                    self.refs[starting_chunk_idx] = starting_chunk;
                    self.refs.insert(starting_chunk_idx + 1, ending_chunk);
                },
            };

            return Some(());
        }

        // If we are not adding an entry, then all of the operations will be done in-place.

        let mut running_length = starting_chunk_idx;

        // Adds the chunks
        if starting_chunk.len() != 0 {
            self.refs[running_length] = starting_chunk;
            running_length += 1;
        }

        if ending_chunk.len() != 0 {
            self.refs[running_length] = ending_chunk;
            running_length += 1;
        }

        // unsafe { ptr::copy::<T>(self.refs[ending_chunk_idx + 1] as *const _, self.refs[running_length] as *mut _, self.refs.len() - ending_chunk_idx - 1) };
        // This might be able to be replaced with a [`ptr::copy`] call however, in many cases this
        // might just be shifting one element at a time where the speedups may be very little.
        // We subtract 1 from the length because of the padded value added earlier.
        for i in (ending_chunk_idx + 1)..(self.refs.len()) {
            self.refs[running_length] = self.refs[i];
            running_length += 1;
        }

        // running_length += self.refs.len() - ending_chunk_idx - 2;
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
    /// # use std::sync::Arc;
    /// // Creating a sevec
    /// let mut sevec: Sevec<u32> = vec![1, 2, 3].into();
    ///
    /// // Creating new data
    /// let slice_data = Arc::new([4, 5, 6]);
    ///
    /// // Adding the data before anything else.
    /// sevec.insert_arc_slice_to_chunk_pos(0, slice_data);
    ///
    /// // Checking the result.
    /// assert_eq!(&sevec.to_string(), "[4, 5, 6, 1, 2, 3]");
    /// ```
    pub fn insert_arc_slice_to_chunk_pos(&mut self, chunk_index: usize, value: Arc<[T]>) -> () {
        // Gets the reference
        let data_inner_ref = ptr::slice_from_raw_parts(value.as_ptr(), value.len());
        // Adds the data.
        self.data.push(value);
        // Adds the reference.
        self.refs.insert(chunk_index, data_inner_ref);
        return ();
    }

    /// Gets the inner reference values.
    pub fn get_inner_refs(&self) -> &[*const [T]] {
        return &self.refs;
    }

    /// Gets the inner reference values mutably.
    /// Adding values that aren't in scope for the lifetime of this object will cause undefined
    /// behavior.
    /// Know what you're doing before using this method.
    pub unsafe fn get_inner_refs_mut(&mut self) -> &mut Vec<*const [T]> {
        return &mut self.refs;
    }

    /// Gets the inner data values.
    pub fn get_inner_data(&self) -> &[Arc<[T]>] {
        return &self.data;
    }

    /// Gets the inner data values mutably.
    /// Removing values that aren't in scope for the lifetime of this object will cause undefined
    /// behavior.
    /// Know what you're doing before using this method.
    pub unsafe fn get_inner_data_mut(&mut self) -> &mut Vec<Arc<[T]>> {
        return &mut self.data;
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

impl <T> From<Arc<[T]>> for Sevec<T> {
    fn from(value: Arc<[T]>) -> Self {
        // Gets the length
        let value_len = value.len();
        let ptr = ptr::slice_from_raw_parts(value.as_ptr(), value_len);
        return Self {
            data: vec![ value, ],
            refs: vec![ ptr,   ],
        };
    }
}

impl <T: Clone + Sized> From<&[T]> for Sevec<T> {
    fn from(value: &[T]) -> Self {
        let mut data = Self::new();
        data.push_slice(value);
        return data;
    }
}

impl <T: Clone + Sized> From<Vec<T>> for Sevec<T> {
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

#[cfg(feature = "serde")]
mod serde_impl {

    use super::*;

    #[derive(Default)]
    pub struct SevecVisitor;

    impl <'de> serde::de::Visitor<'de> for SevecVisitor {
        type Value = Sevec<u8>;
        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            return formatter.write_str("Some Bytes.");
        }
        fn visit_bytes<E>(self, v: &[u8]) -> Result<Self::Value, E>
            where
                E: serde::de::Error, {
            let res = Self::Value::from(v);
            return Ok(res);
        }
        fn visit_byte_buf<E>(self, v: Vec<u8>) -> Result<Self::Value, E>
            where
                E: serde::de::Error, {
            return Ok(v.into());
        }

        fn visit_borrowed_bytes<E>(self, v: &'de [u8]) -> Result<Self::Value, E>
            where
                E: serde::de::Error, {
            return Ok(Self::Value::from(v));
        }

        fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: serde::de::SeqAccess<'de>, {

            // As it should be but implementations don't like me.
            // let len = seq.size_hint().ok_or(serde::de::Error::custom("Failed to get `size_hint` for `Sevec<T>` visitor!"))?;

            let res = match seq.size_hint() {
                // Better path with pre-allocation with size-hint.
                Some(len) => {
                    // Gets data
                    let mut data = {
                        let data = Arc::<[u8]>::new_uninit_slice(len);
                        unsafe { data.assume_init() }
                    };
                    let data_mut = Arc::get_mut(&mut data).unwrap();
                    // Gets running index
                    let mut count = 0;
                    // Writes the data
                    while let Ok(Some(v)) = seq.next_element() {

                        // This is a bad sign!
                        // We just break in this case but this does mean that the size hint
                        // lied to us.
                        if count >= len {
                            return Err(serde::de::Error::custom("size_hint doesn't match actual size (size_hint undershot), a full array cannot be created."));
                        }
                        // Updates the value
                        data_mut[count] = v;
                        count += 1;
                    }
                    // let mut value = Sevec::new();
                    // value.push_arc_slice(data);
                    // return value;
                    Self::Value::from(data)
                },

                // Worse path if we don't know the length.
                None => {
                    let mut buf = Vec::new();
                    while let Ok(Some(v)) = seq.next_element() {
                        buf.push(v);
                    }

                    Self::Value::from(buf)
                },

            };

            return Ok(res);

        }


    }

    impl <'de> serde::Deserialize<'de> for Sevec<u8> {
        fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
            where
                D: serde::Deserializer<'de> {
            return deserializer.deserialize_bytes(SevecVisitor);
        }

    }

    impl serde::Serialize for Sevec<u8> {
        fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
            where
                S: serde::Serializer {
            let data: Vec<_> = self.into();
            return serializer.serialize_bytes(&data);
        }
    }

    #[test]
    fn test_serde() {

        let data: Sevec<_> = vec![1, 2, 3, 4].into();

        // Serializes
        let string = serde_json::to_string(&data).unwrap();
        // Deserializes
        let deser: Sevec<u8> = serde_json::from_str(&string).unwrap();

        // Is the same as the original
        assert_eq!(deser.to_string(), data.to_string());

    }

}

// SAFETY: Raw pointers reference locally owned `Arc<T>` data.
// Sound if `T: Send + Sync`.
// This impl is copied from Arc<T>'s impl of both of these traits.
unsafe impl <T: Send + Sync> Send for Sevec<T> {}
unsafe impl <T: Send + Sync> Sync for Sevec<T> {}

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
        v.insert_arc_slice_to_chunk_pos(0, vec!["Hello", "H"].into());
        v.remove_range(0..2).unwrap();
        assert_eq!(format!("{:?}", vec), format!("{:?}", v));
    }

    #[test]
    fn test_insert_slice() {

        let mut sevec = Sevec::new();

        assert!(sevec.insert_slice(1, &[1, 2, 3]).is_none());
        assert_eq!(sevec.len(), 0);
        assert_eq!(sevec.refs.len(), 0);

        sevec.insert_slice(0, &[1, 2, 3]).unwrap();
        assert_eq!(&sevec.to_string(), "[1, 2, 3]");

        // Testing adding a value to the end.
        sevec.insert_slice(3, &[4, 5, 6]).unwrap();
        assert_eq!(&sevec.to_string(), "[1, 2, 3, 4, 5, 6]");

        // Testing adding a value inside a slice.
        sevec.insert_slice(1, &[7, 8, 9]).unwrap();
        assert_eq!(&sevec.to_string(), "[1, 7, 8, 9, 2, 3, 4, 5, 6]");

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

        // Testing adding a value to the end.
        unsafe { sevec.insert_raw_slice(6, data_ptr) }.unwrap();
        assert_eq!(&sevec.to_string(), "[4, 5, 6, 4, 5, 6, 4, 5, 6]");

        // Testing adding a value inside a slice.
        unsafe { sevec.insert_raw_slice(1, data_ptr) }.unwrap();
        assert_eq!(&sevec.to_string(), "[4, 4, 5, 6, 5, 6, 4, 5, 6, 4, 5, 6]");

        // Testing adding a value at the start.
        unsafe { sevec.insert_raw_slice(0, data_ptr) }.unwrap();
        assert_eq!(&sevec.to_string(), "[4, 5, 6, 4, 4, 5, 6, 5, 6, 4, 5, 6, 4, 5, 6]");

    }

    #[test]
    fn test_remove_basic() {

        let mut sevec = Sevec::new();
        sevec.push(1);
        sevec.push(2);
        sevec.push(3);
        sevec.push(4);

        assert_eq!(&sevec.to_string(), "[1, 2, 3, 4]");
        sevec.remove_range(1..=2).unwrap();

        assert_eq!(sevec.len(), 2);
        assert_eq!(&sevec.to_string(), "[1, 4]");

        sevec.remove(1).unwrap();
        assert_eq!(&sevec.to_string(), "[1]");

        sevec.insert_slice(1, &[5, 6, 7, 8]);
        assert_eq!(&sevec.to_string(), "[1, 5, 6, 7, 8]");

        sevec.remove(2);
        assert_eq!(&sevec.to_string(), "[1, 5, 7, 8]");

        sevec.insert_slice(1, &[9, 10, 11]);
        assert_eq!(&sevec.to_string(), "[1, 9, 10, 11, 5, 7, 8]");

    }

    #[test]
    fn test_remove_edge_case_1() {

        let mut sevec = Sevec::new();
        sevec.push(1);
        assert_eq!(&sevec.to_string(), "[1]");

        sevec.insert_slice(1, &[2, 3, 4]);
        assert_eq!(&sevec.to_string(), "[1, 2, 3, 4]");

        sevec.push(5);
        assert_eq!(&sevec.to_string(), "[1, 2, 3, 4, 5]");

        sevec.remove(2);
        assert_eq!(&sevec.to_string(), "[1, 2, 4, 5]");

    }

    #[test]
    fn test_remove_edge_case_2() {

        let mut sevec = Sevec::new();

        sevec.insert_slice(0, &[1]);
        sevec.insert_slice(1, &[5, 6, 7]);
        sevec.insert_slice(1, &[2, 3, 4]);

        assert_eq!(&sevec.to_string(), "[1, 2, 3, 4, 5, 6, 7]");

        sevec.remove(3).unwrap();
        assert_eq!(&sevec.to_string(), "[1, 2, 3, 5, 6, 7]");

    }

    #[test]
    fn test_remove_edge_case_3() {

        let mut sevec = Sevec::new();

        sevec.insert_slice(0, &[1, 2, 3, 4]).unwrap();
        assert_eq!(&sevec.to_string(), "[1, 2, 3, 4]");

        sevec.insert_slice(3, &[0]).unwrap();
        assert_eq!(&sevec.to_string(), "[1, 2, 3, 0, 4]");

        sevec.remove(1).unwrap();
        assert_eq!(&sevec.to_string(), "[1, 3, 0, 4]");
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
        data.push_arc_slice(vec![5, 2, 1].into());
        data.push(1);
        data.push(2);
        data.push(3);
        data.push(4);
        data.push_arc_slice(vec![1, 2, 3].into());

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
        data.push_arc_slice(vec![1, 2, 3].into());

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
    fn test_remove_and_copy_from_end() {

        let mut data_1: Sevec<u8> = Sevec::new();
        data_1.push(1);
        data_1.push(2);
        data_1.push(3);
        data_1.push(4);

        assert!(data_1.remove_and_copy_slice_from_end(5).is_none());

        let res_1 = data_1.remove_and_copy_slice_from_end(2).unwrap();

        assert_eq!(data_1.to_string(), "[1, 2]");
        assert_eq!(format!("{:?}", res_1), "[3, 4]");

        let res_2 = data_1.remove_and_copy_slice_from_end(0).unwrap();
        assert_eq!(data_1.to_string(), "[1, 2]");
        assert_eq!(format!("{:?}", res_2), "[]");

        let res_3 = data_1.remove_and_copy_slice_from_end(1).unwrap();
        assert_eq!(data_1.to_string(), "[1]");
        assert_eq!(format!("{:?}", res_3), "[2]");

        let mut data_2: Sevec<u8> = Sevec::new();
        data_2.push_slice(&[1, 2]);
        data_2.push(3);


        let res_4 = data_2.remove_and_copy_slice_from_end(2).unwrap();
        assert_eq!(data_1.to_string(), "[1]");
        assert_eq!(format!("{:?}", res_4), "[2, 3]");

    }

    #[test]
    fn test_get_refs() {

        let get_elem = |v: &Sevec<usize>, elem: usize| -> Option<usize> {
            let flattened_data = v.get_refs().into_iter().map(|v| *v).flatten().collect::<Vec<_>>();
            let elem = flattened_data.get(elem);
            match elem {
                Some(&&v) => Some(v),
                None => None,
            }
        };

        let mut data = Sevec::new();

        // With Data Checks
        data.push(0);
        data.push(1);
        data.push(2);

        assert_eq!(get_elem(&data, 0), Some(0));
        assert_eq!(get_elem(&data, 1), Some(1));
        assert_eq!(get_elem(&data, 2), Some(2));
        assert_eq!(get_elem(&data, 3), None);

        // Pushes new data
        data.push_arc_slice(vec![1, 2, 3].into());

        // Previous Values
        assert_eq!(get_elem(&data, 0), Some(0));
        assert_eq!(get_elem(&data, 1), Some(1));
        assert_eq!(get_elem(&data, 2), Some(2));

        // Extended values
        assert_eq!(get_elem(&data, 3), Some(1));
        assert_eq!(get_elem(&data, 4), Some(2));
        assert_eq!(get_elem(&data, 5), Some(3));
        assert_eq!(get_elem(&data, 6), None);

        // Data removed check
        data.remove(1).unwrap();

        assert_eq!(get_elem(&data, 0), Some(0));
        assert_eq!(get_elem(&data, 1), Some(2));
        assert_eq!(get_elem(&data, 5), None); // After removal the end would have moved.

        data.remove_range(..).unwrap();
        assert_eq!(data.len(), 0);

    }

    #[test]
    fn test_insert_slice_edge_cases() {

        let mut sevec = Sevec::new();
        sevec.insert_slice(0, &[1, 2, 3]);
        sevec.insert_slice(1, &[]);

        assert_eq!(&sevec.to_string(), "[1, 2, 3]");
    }

    #[test]
    fn test_simple_remove() {

        let mut sevec: Sevec<_> = vec![1, 2, 3].into();
        assert_eq!(&sevec.to_string(), "[1, 2, 3]");

        sevec.remove(1).unwrap();
        // panic!("{:?}", sevec.get_refs());
        assert_eq!(&sevec.to_string(), "[1, 3]");


        sevec.remove(1).unwrap();
        assert_eq!(&sevec.to_string(), "[1]");

    }

}
