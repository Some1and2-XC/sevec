use sevec::*;

pub fn main() {

    let mut v = Sevec::new();
    // v.push("Hello There!");

    // v.push("Hello There!");
    // v.push("Hello There!");
    // v.push("Hello There!");
    v.push("BBB");

    v.push("\0");
    // let long_string = (0..(16384 * 16384)).into_iter().map(|_v| 'A').collect::<String>();
    let long_string = (0..16384).into_iter().map(|_v| 'A').collect::<String>();
    v.push(&long_string);

    // println!("Here is the array: {:?}", v);

    let sevec_size = v.size_estimation();
    let vec_size = size_of::<Vec<()>>();

    println!("\tCopy Statistics:");

    println!("Here is the estimated size of the array: {}", sevec_size);
    // println!("Here is the estimated size of a normal array: {}", vec_size);
    // println!("The array is {:.2}x more than a normal array!", (sevec_size as f64 / vec_size as f64));

    let mut bytes_total = 0;
    let mut i = 0;

    while let Some(chunk) = v.get_chunk(i) {
        bytes_total += chunk.iter().map(|v| v.len()).sum::<usize>();
        i += 1;
    }

    println!("The inner data however is {} bytes which works out to {} bytes to clone", bytes_total, vec_size + bytes_total);
    println!("This works out to {:.2}x more than a normal array (which is the amount to copy)!", ((vec_size + bytes_total) as f64) / sevec_size as f64);

}
