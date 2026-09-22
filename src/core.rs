pub(crate) mod grid;
pub(crate) mod irradiance;
pub(crate) mod solar_position;
pub(crate) mod time;
pub(crate) mod utils;

use rayon::{ThreadPool, ThreadPoolBuilder};
use std::collections::HashMap;
use std::sync::{Arc, OnceLock, RwLock};

use crate::core::grid::SolarError;

fn get_or_create_thread_pool(num_threads: usize) -> Result<Arc<ThreadPool>, SolarError> {
    static POOL_CACHE: OnceLock<RwLock<HashMap<usize, Arc<ThreadPool>>>> = OnceLock::new();

    let cache = POOL_CACHE.get_or_init(|| RwLock::new(HashMap::new()));

    {
        let reader = cache
            .read()
            .map_err(|e| SolarError::ThreadPool(e.to_string()))?;
        if let Some(pool) = reader.get(&num_threads) {
            return Ok(Arc::clone(pool));
        }
    }

    let mut writer = cache
        .write()
        .map_err(|e| SolarError::ThreadPool(e.to_string()))?;
    {
        if let Some(pool) = writer.get(&num_threads) {
            return Ok(Arc::clone(pool));
        }
    }

    let pool = Arc::new(
        ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .map_err(|e| SolarError::ThreadPool(e.to_string()))?,
    );

    writer.insert(num_threads, Arc::clone(&pool));
    Ok(pool)
}

pub(crate) fn with_thread_pool<F, R>(num_threads: usize, op: F) -> Result<R, SolarError>
where
    F: FnOnce() -> R + Send,
    R: Send,
{
    let pool = get_or_create_thread_pool(num_threads)?;
    Ok(pool.install(op))
}
