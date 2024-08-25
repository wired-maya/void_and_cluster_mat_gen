use core::f64;
use std::{env, io::{stdout, Cursor, Write}};
use image::{ImageError, Rgb, RgbImage};

const SIGMA: f64 = 1.5;

fn main() -> Result<(), ImageError> {
    let args: Vec<String> = env::args().collect();

    let m: usize = *(&args[1].parse().expect("Arguments 1 and 2 (M, N) should be positive integers"));
    let n: usize = *(&args[2].parse().expect("Arguments 1 and 2 (M, N) should be positive integers"));

    // Generate input pattern, an MxN matrix where <= 50% of the pixels are 1, the rest being 0
    // Do this based on random noise
    let input_pattern: Vec<usize> = gen_input_pattern(m, n);

    // Evenly distribute the pixels within the input pattern
    let initial_bin_pattern = gen_initial_binary_pattern(&input_pattern, m, n);

    // Convert to dither pattern
    let dither_pattern: Vec<usize> = gen_dither_array(&initial_bin_pattern, m, n);

    // Convert to image buffer
    let mut img = RgbImage::new(m as u32, n as u32);

    for i in 0..n {
        for j in 0..m {
            let original_val: f32 = dither_pattern[i * m + j] as f32;
            let new_val: f32 = (original_val / (m * n) as f32) * 255.0;
            let b: u8 = new_val as u8;

            img.put_pixel(j as u32, i as u32, Rgb([b, b, b]));
        }
    }

    // Write bytes to console
    let mut bytes: Vec<u8> = Vec::new();
    img.write_to(&mut Cursor::new(&mut bytes), image::ImageFormat::Png)?;
    stdout().write_all(&bytes)?;

    Ok(())
}

fn _print_mat_bin(mat: &Vec<usize>, m: usize, n: usize) {
    for i in 0..n {
        for j in 0..m {
            if mat[i * m + j] == 1 {
                print!("#");
            } else {
                print!(" ");
            }
            
        }

        print!("\n\n");
    }
}

fn _print_mat_f64(mat: &Vec<f64>, m: usize, n: usize) {
    for i in 0..n {
        for j in 0..m {
            print!("{:.2} ", mat[i * m + j]);
        }

        print!("\n\n");
    }
}

fn _print_mat_dither(mat: &Vec<usize>, m: usize, n: usize) {
    for i in 0..n {
        for j in 0..m {
            print!("{:0>3} ", mat[i * m + j]);
        }

        print!("\n\n");
    }
}

// Generate random white noise pattern
// In this case since black pixels (1) are <= 50%, they are the minority pixels
fn gen_input_pattern(m: usize, n: usize) -> Vec<usize> {
    let mut vec: Vec<usize> = Vec::new();
    let vec_len: usize = m * n;
    let mut black_count: usize = 0;

    // Generate random values in vec
    for _ in 0..vec_len {
        let mut val = 0;
        if rand::random() {
            val = 1;
            black_count += 1;
        }

        vec.push(val);
    }

    // If 1s are > 50%, invert colours
    if black_count > vec_len / 2 {
        for i in 0..vec_len as usize {
            if vec[i] == 1 {
                vec[i] = 0;
            } else {
                vec[i] = 1;
            }
        }
    }

    vec
}

// Generate weight map that represents how clustered each pixel is
fn gen_weight_map(mat: &Vec<usize>, m: usize, n: usize) -> Vec<f64> {
    let mut weight_map: Vec<f64> = Vec::new();

    for i in 0..n {
        for j in 0..m {
            let mut weight: f64 = 0.0;

            // Calculate contribution of weight of each pixel
            // Treat current pixel as centre, so wrap around coords
            for k in -(n as isize / 2)..(n as isize / 2) {
                for l in -(m as isize / 2)..(m as isize / 2) {
                    let abs_coords_x: usize = (((m + j) as isize - l) % m as isize) as usize;
                    let abs_coords_y: usize = (((n + i ) as isize - k) % n as isize) as usize;
                    
                    let d_sq: f64 = (k * k + l * l) as f64;

                    let weight_contrib: f64 = f64::consts::E.powf(- d_sq / (2.0 * SIGMA * SIGMA));

                    weight += weight_contrib * mat[abs_coords_y * n + abs_coords_x] as f64;
                }
            }

            weight_map.push(weight);
        }
    }

    weight_map
}

fn find_largest_void(weight_map: &Vec<f64>, bin_map: &Vec<usize>, m: usize, n: usize) -> (usize, usize) {
    let mut weight_least: f64 = f64::MAX;
    let mut weight_least_x: usize = 0;
    let mut weight_least_y: usize = 0;

    for i in 0..n {
        for j in 0..m {
            if bin_map[i * m + j] == 0 { // Only 0s are under considerations
                let weight = weight_map[i * m + j];

                if weight < weight_least {
                    weight_least = weight;
                    weight_least_x = j;
                    weight_least_y = i;
                }
            }
        }
    }

    (weight_least_x, weight_least_y)
}

fn find_tightest_cluster(weight_map: &Vec<f64>, bin_map: &Vec<usize>, m: usize, n: usize) -> (usize, usize) {
    let mut weight_greatest: f64 = f64::MIN;
    let mut weight_greatest_x: usize = 0;
    let mut weight_greatest_y: usize = 0;

    for i in 0..n {
        for j in 0..m {
            if bin_map[i * m + j] == 1 { // Only 1s are under considerations
                let weight = weight_map[i * m + j];

                if weight > weight_greatest {
                    weight_greatest = weight;
                    weight_greatest_x = j;
                    weight_greatest_y = i;
                }
            }
        }
    }

    (weight_greatest_x, weight_greatest_y)
}

fn gen_initial_binary_pattern(mat: &Vec<usize>, m: usize, n: usize) -> Vec<usize> {
    // 1. Find location of tightest cluster (all 1s are candidates)
    // 2. Remove 1 with the tightest cluster
    // 3. Find location of largest void (all 0s are candidates)
    // 4. Did removing the 1 create the largest void?
    // 5a. Yes: insert a 1 in the largest void and jump to step 1
    // 5b. No: Restore 1 that was just removed and return

    let mut mat: Vec<usize> = mat.clone();

    loop {
        let tc_coords: (usize, usize) = find_tightest_cluster(&gen_weight_map(&mat, m, n), &mat, m, n);
        
        mat[tc_coords.1 * n + tc_coords.0] = 0;

        let lv_coords: (usize, usize) = find_largest_void(&gen_weight_map(&mat, m, n), &mat, m, n);

        mat[lv_coords.1 * n + lv_coords.0] = 1;

        if tc_coords == lv_coords { break; }
    }

    mat
}

fn gen_dither_array(initial_binary_pattern: &Vec<usize>, m: usize, n: usize) -> Vec<usize> {
    // Phase 1 Start
    //      prototype binary pattern = initial binary pattern
    //      ones = 1 count in prototype binary pattern
    //      rank = ones - 1
    // Phase 1 Loop
    //      If rank is < 0, go to phase 2
    //      Find tightest cluster
    //      Change pixel to 0
    //      Enter current rank value into corresponding location in dither array
    //      rank -= 1
    // Phase 2 Start
    //      prototype binary pattern = initial binary pattern
    //      rank = ones
    // Phase 2 Loop
    //      If rank is >= MN / 2, go to phase 3
    //      Find largest void
    //      Change pixel to 1
    //      Enter current rank value into corresponding location in dither array
    //      rank += 1;
    // Phase 3 Start
    //      invert prototype binary pattern
    // Phase 3 Loop
    //      If rank >= MN, return
    //      Find tightest cluster
    //      Change pixel to 0
    //      Enter current rank value into corresponding location in dither array
    //      rank += 1

    let mut dither_array: Vec<usize> = vec![0; m * n];

    // Phase 1 Start
    let mut prototype_bin_pat: Vec<usize> = initial_binary_pattern.clone();
    let mut ones: isize = 0;

    for i in 0..m * n {
        if prototype_bin_pat[i] == 1 { ones +=1; }
    }

    let mut rank: isize = ones - 1;

    // Phase 1 Loop
    loop {
        if rank < 0 { break; }

        let tc_coords: (usize, usize) =
            find_tightest_cluster(&gen_weight_map(&prototype_bin_pat, m, n), &prototype_bin_pat, m, n);
        
        prototype_bin_pat[tc_coords.1 * n + tc_coords.0] = 0;
        dither_array[tc_coords.1 * n + tc_coords.0] = rank as usize; // Can't be -1 so casting always works

        rank -= 1;
    }

    // Phase 2 Start
    prototype_bin_pat = initial_binary_pattern.clone();
    rank = ones;

    // Phase 2 Loop
    loop {
        if rank >= ((m * n) / 2) as isize { break; }

        let lv_coords: (usize, usize) =
            find_largest_void(&gen_weight_map(&prototype_bin_pat, m, n), &prototype_bin_pat, m, n);
        
        prototype_bin_pat[lv_coords.1 * n + lv_coords.0] = 1;
        dither_array[lv_coords.1 * n + lv_coords.0] = rank as usize;
        rank += 1;
    }

    // Phase 3 Start
    for i in 0..(m * n) as usize {
        if prototype_bin_pat[i] == 1 {
            prototype_bin_pat[i] = 0;
        } else {
            prototype_bin_pat[i] = 1;
        }
    }

    // Phase 3 Loop
    loop {
        if rank >= (m * n) as isize { break; }

        let tc_coords: (usize, usize) =
            find_tightest_cluster(&gen_weight_map(&prototype_bin_pat, m, n), &prototype_bin_pat, m, n);
        
        prototype_bin_pat[tc_coords.1 * n + tc_coords.0] = 0;
        dither_array[tc_coords.1 * n + tc_coords.0] = rank as usize; // Can't be -1 so casting always works

        rank += 1;
    }

    dither_array
}