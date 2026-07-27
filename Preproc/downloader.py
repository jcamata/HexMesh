import os
import sys
import time
import gzip
import shutil

def _format_time(seconds: float) -> str:
    if seconds < 0 or seconds == float('inf'):
        return "--:--"
    mins, secs = divmod(int(seconds), 60)
    hrs, mins = divmod(mins, 60)
    if hrs > 0:
        return f"{hrs:02d}:{mins:02d}:{secs:02d}"
    return f"{mins:02d}:{secs:02d}"

def download_file_if_missing(url: str, local_path: str, config=None) -> bool:
    if os.path.exists(local_path):
        print(f"[CACHE] File '{local_path}' already exists. Skipping download.")
        return True

    if not url:
        print(f"[WARNING] No download URL provided for '{local_path}'. Skipping.")
        return False

    if config:
        lat_first = max(config.lat_min, config.lat_max)
        lat_last = min(config.lat_min, config.lat_max)
        url = url.format(
            lat_min=config.lat_min,
            lat_max=config.lat_max,
            lon_min=config.lon_min,
            lon_max=config.lon_max,
            lat_first=lat_first,
            lat_last=lat_last,
            north=config.lat_max,
            south=config.lat_min,
            west=config.lon_min,
            east=config.lon_max
        )

    filename = os.path.basename(local_path)
    print(f"\n[DOWNLOAD] File not found locally: {local_path}")
    print(f"[DOWNLOAD] Connecting to: {url}")

    is_gz_download = url.endswith('.gz') and not local_path.endswith('.gz')
    # A .zip URL for a file inside it (e.g. GSHHG shapefiles): download the archive next to
    # the requested file and unpack it there. local_path names a file within the archive tree.
    is_zip_download = url.endswith('.zip') and not local_path.endswith('.zip')
    zip_root = None
    if is_zip_download:
        zip_root = os.path.dirname(local_path.rstrip(os.sep)).split(os.sep)[0] or '.'
        target_write_path = zip_root + '.zip'
    else:
        target_write_path = local_path + '.gz' if is_gz_download else local_path

    try:
        import requests
        try:
            from tqdm import tqdm
            HAS_TQDM = True
        except ImportError:
            HAS_TQDM = False

        start_time = time.time()
        headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7)'}
        response = requests.get(url, stream=True, headers=headers, timeout=90)
        response.raise_for_status()

        total_size = int(response.headers.get('content-length', 0))
        block_size = 1024 * 64

        if total_size > 0:
            print(f"[DOWNLOAD] Total File Size: {total_size / (1024 * 1024):.2f} MB")

        with open(target_write_path, 'wb') as f:
            if HAS_TQDM:
                with tqdm(
                    desc=f"Download {filename}",
                    total=total_size if total_size > 0 else None,
                    unit='B',
                    unit_scale=True,
                    unit_divisor=1024,
                    ncols=90,
                    ascii=True,
                    bar_format="{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]" if total_size > 0 else "{desc}: {n_fmt} [{elapsed}, {rate_fmt}]"
                ) as bar:
                    for data in response.iter_content(block_size):
                        size = f.write(data)
                        bar.update(size)
            else:
                downloaded = 0
                last_print_time = time.time()

                for data in response.iter_content(block_size):
                    size = f.write(data)
                    downloaded += size
                    now = time.time()

                    if now - last_print_time >= 0.25 or (total_size > 0 and downloaded >= total_size):
                        elapsed = now - start_time
                        speed = downloaded / elapsed if elapsed > 0 else 0
                        speed_mb = speed / (1024 * 1024)
                        dl_mb = downloaded / (1024 * 1024)

                        if total_size > 0:
                            percent = (downloaded / total_size) * 100
                            tot_mb = total_size / (1024 * 1024)
                            eta_sec = (total_size - downloaded) / speed if speed > 0 else 0
                            bar_len = 25
                            filled = int(bar_len * downloaded // total_size)
                            bar = '=' * filled + '>' + '.' * (bar_len - filled - 1)
                            sys.stdout.write(
                                f"\r[DOWNLOAD] {filename}: [{bar}] {percent:5.1f}% | {dl_mb:.1f}/{tot_mb:.1f} MB | {speed_mb:.2f} MB/s | ETA: {_format_time(eta_sec)}"
                            )
                        else:
                            sys.stdout.write(
                                f"\r[DOWNLOAD] {filename}: {dl_mb:.1f} MB downloaded | {speed_mb:.2f} MB/s | Elapsed: {_format_time(elapsed)}"
                            )
                        sys.stdout.flush()
                        last_print_time = now

                print()

        if is_gz_download:
            print(f"[DOWNLOAD] Decompressing gzip archive '{target_write_path}' to '{local_path}'...")
            with gzip.open(target_write_path, 'rb') as f_in:
                with open(local_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(target_write_path)

        if is_zip_download:
            import zipfile
            print(f"[DOWNLOAD] Extracting '{target_write_path}' into '{zip_root}/'...")
            with zipfile.ZipFile(target_write_path) as z:
                z.extractall(zip_root)
            os.remove(target_write_path)
            if not os.path.exists(local_path):
                print(f"[ERROR] '{local_path}' not found in the archive.")
                return False

        print(f"[DOWNLOAD] Complete: '{local_path}'\n")
        return True
    except Exception as e:
        print(f"\n[ERROR] Download failed for '{local_path}': {e}")
        if os.path.exists(target_write_path):
            os.remove(target_write_path)
        if os.path.exists(local_path):
            os.remove(local_path)
        return False
