// Simple program for finding the Weibull distribution parameters
// k shape factor and c scale factor.
// The file with the example wind measurement data set comes
// from the measurement mast US Virgin Islands St. Thomas Bovoni and
// was downloaded from the site 
// https://midcdmz.nrel.gov/apps/sitehome.pl?site=USVILONA.
// 
// Roberts, O.; Andreas, A.; (1997). United States Virgin Islands:
// St. Thomas & St. Croix (Data); NREL Report No. DA-5500-64451.
// http://dx.doi.org/10.7799/1183464 
// https://midcdmz.nrel.gov/
//
// https://midcdmz.nrel.gov/
//
// Sorting function from 
// https://fortran-lang.discourse.group/t/modern-fortran-sample-code/2019/4
//

use std::env;
use std::error::Error;
use std::ffi::OsString;
use std::process;

use serde::Deserialize;

#[derive(Deserialize)]
struct WindData {
    // wind speed values
    wind_speeds: f64,
}

const TAYLOR_COEFFICIENTS: [f64; 29] = [
    -0.00000000000000000023,  0.00000000000000000141,  0.00000000000000000119,
    -0.00000000000000011813,  0.00000000000000122678, -0.00000000000000534812,
    -0.00000000000002058326,  0.00000000000051003703, -0.00000000000369680562,
     0.00000000000778226344,  0.00000000010434267117, -0.00000000118127457049,
     0.00000000500200764447,  0.00000000611609510448, -0.00000020563384169776,
     0.00000113302723198170, -0.00000125049348214267, -0.00002013485478078824,
     0.00012805028238811619, -0.00021524167411495097, -0.00116516759185906511,
     0.00721894324666309954, -0.00962197152787697356, -0.04219773455554433675,
     0.16653861138229148950, -0.04200263503409523553, -0.65587807152025388108,
     0.57721566490153286061,  1.00000000000000000000,
];

const INITIAL_SUM: f64 = 0.00000000000000000002;

fn main() {
    if let Err(err) = find_weibull() {
        println!("{}", err);
        process::exit(1);
    }
}

fn get_first_arg() -> Result<OsString, Box<dyn Error>> {
    match env::args_os().nth(1) {
        None => Err(From::from("Stop! Usage cargo run <filename>")),
        Some(file_path) => Ok(file_path),
    }
}

// https://codereview.stackexchange.com/questions/237790/quick-sort-algorithm-in-rust
fn quick_sort<T: PartialOrd>(a: &mut [T]) {
    let n = a.len();
    if n > 1 {
        let mut j = 0;
        for i in 0..n {
            if &a[i] < &a[n - 1] {
                a.swap(i, j);
                j += 1;
            }
        }
        a.swap(j, n - 1); // pivot
        quick_sort(&mut a[0..j]);
        quick_sort(&mut a[j + 1..n]);
    }
}

fn median(x: &Vec<f64>) -> f64 {
    
    let mut xs = x.clone();
    
    quick_sort(&mut xs);
    
    let n: usize = xs.len();

    if n % 2 == 0 {
        let i  = n/2;
        let i2 = i + 1;
        let res: f64 = (xs[i] + xs[i2]) / 2.0;
        return res;
    } else {
        let i = (n+1) / 2;
        let res: f64 = xs[i];
        return res;
    }
}

fn gamma(x: f64) -> f64 {
    TAYLOR_COEFFICIENTS.iter().fold(INITIAL_SUM, |sum, coefficient| {
        sum * (x - 1.0) + coefficient
    }).recip()
}

fn k_estimator(x: &Vec<f64>, kin: &f64) -> f64 {
    
    let n = x.len() as f64;
    
    let xc = x.clone();
    let mut x3 = Vec::new();
    
    // iterate immutably
    // elements are immutable pointers
    for item in &xc {
        let item: &f64 = item;
        x3.push(item.powf(3.0));
    }
    let sum1: f64 = x3.iter().sum::<f64>() / n;
    let sum2: f64 = xc.iter().sum::<f64>() / n;
    let sum2: f64 = sum2 * sum2 * sum2;
    let gamma1: f64 = gamma(1.0+1.0/kin);
    let gamma1: f64 = gamma1 * gamma1 * gamma1;
    let gamma2: f64 = gamma(1.0+3.0/kin);
    let res: f64 = sum1 / sum2;
    let res: f64 = res * gamma1 - gamma2;
    return res
}

fn bisection(x: &Vec<f64>, ikmin: f64, ikmax: f64, eps: f64, iter: u32) -> f64 {
    
    let xs = x.clone();
    
    // initial values
    let mut kmin = ikmin.clone();
    let mut kmax = ikmax.clone();
    let mut fkmin: f64 = k_estimator(&xs, &kmin);
    let mut fkmax: f64 = k_estimator(&xs, &kmax);
    
    if fkmin * fkmax > 0.0 {
        println!("Error: Both estimated k values are greater than zero!");
        return 0.0
    }
    
    for _j in 1..iter {
        let k: f64 = (kmin + kmax) / 2.0;
        let fk: f64 = k_estimator(&xs, &k);
        let fkk: f64 = (fkmax - fkmin) / (kmax - kmin);
        
        let fk_fkk: f64 = fk/fkk;
        if fk_fkk.abs() - eps > 0.0 {
            if fk*fkmin < 0.0 {
                kmax = k;
                fkmax = fk;
            } else {
                if fk * fkmin == 0.0 {
                    return k
                }
                kmin = k;
                fkmin = fk;
            }
        } else {
            return k
        }
    }
    return 0.0
}
    
    
fn find_weibull() -> Result<(), Box<dyn Error>> {
    
    // vector for storing wind speed data
    let mut ws = Vec::new();
    
    // get the first argument 
    let file_path = get_first_arg()?;
    
    let mut rdr = csv::Reader::from_path(file_path)?;
    
    for result in rdr.deserialize() {
        let record = WindData {
            wind_speeds: result?,
        };
        ws.push(record.wind_speeds);
    }
    
    // mean wind speed
    let n = ws.len() as f64;
    let ws_mean: f64 = ws.iter().sum::<f64>() / n;
    
    // median wind speed
    let ws_median: f64 = median(&ws);
    
    // range for searching k and c
    let kmin: f64 = 1.0;
    let kmax: f64 = 8.0;
    // accuracy
    let eps: f64 = 0.000001;
    // number of iterations
    let niter: u32 = 50;
    
    // find shape factor
    let k: f64 = bisection(&ws, kmin, kmax, eps, niter);
    
    // calculate c scale factor
    let c: f64 = 0.586 + 0.433/k;
    let c: f64 = ws_mean * c.powf(-1.0/k);
    
    println!();
    println!("Found Weibull distribution parameters:");
    println!();
    println!("shape factor k: {:.2}", k);
    println!("scale factor c: {:.2}", c);
    println!();
    println!("Mean wind speed: {:.2} m/s", ws_mean);
    println!("Median wind speed: {:.2} m/s", ws_median);
    println!();

    Ok(())
}
