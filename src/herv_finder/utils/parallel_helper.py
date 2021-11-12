"""Helper for multiprocessing"""
import gc
import multiprocessing
import threading
import time
from typing import Callable

import tqdm


class _DumbTQDM:
    def __init__(self, *args, **kwargs):
        """"""
        pass

    def close(self, *args, **kwargs):
        """"""
        pass

    def update(self, *args, **kwargs):
        """"""
        pass


class MultiThreadingJobQueue(threading.Thread):
    def __init__(self,
                 pool_name: str = "Unnamed pool",
                 pool_size: int = multiprocessing.cpu_count(),
                 with_tqdm: bool = False,
                 return_function: Callable = None):
        super().__init__()
        self.pool_size = pool_size
        self._job_queue = []
        self._is_terminated = False
        self._max_queue_len = 0
        self.with_tqdm = with_tqdm
        self.pool_name = pool_name
        self.return_function = return_function
        self.return_value = []

    def run(self):
        active_thread = []
        if self.with_tqdm:
            pbar = tqdm.tqdm(desc=self.pool_name, total=self._max_queue_len)
        else:
            pbar = _DumbTQDM()
        while len(self._job_queue) > 0 and not self._is_terminated:
            while len(self._job_queue) > 0 and len(active_thread) < self.pool_size:
                active_thread.append(self._job_queue.pop(0))
                active_thread[-1].start()
            for thread in active_thread:
                if not thread.is_alive():
                    if self.return_function is not None:
                        self.return_value.append(self.return_function(thread))
                    thread.join()
                    active_thread.remove(thread)
                    pbar.update(1)
            time.sleep(0.1)
        while len(active_thread) > 0 and not self._is_terminated:
            for thread in active_thread:
                if not thread.is_alive():
                    if self.return_function is not None:
                        self.return_value.append(self.return_function(thread))
                    thread.join()
                    active_thread.remove(thread)
                    pbar.update(1)
            time.sleep(0.1)
        pbar.close()
        self._is_terminated = True

    def __len__(self):
        return len(self._job_queue)

    def stop(self):
        self._is_terminated = True

    def add(self, mt_instance: threading.Thread):
        self._job_queue.append(mt_instance)
        self._max_queue_len += 1

    def append(self, *args, **kwargs):
        self.add(*args, **kwargs)

    def close(self):
        self.join()


class MultiProcessingJobQueue(threading.Thread):
    def __init__(self,
                 pool_name: str = "Unnamed pool",
                 pool_size: int = multiprocessing.cpu_count(),
                 with_tqdm: bool = False,
                 return_function:Callable = None):
        super().__init__()
        self.pool_size = pool_size
        self._job_queue = []
        self._is_terminated = False
        self._max_queue_len = 0
        self.with_tqdm = with_tqdm
        self.pool_name = pool_name

    def run(self):
        def _scan_through_process():
            """
            Scan through all processes and terminate the exited process.
            """
            for process in active_processes:
                if process.exitcode is not None:
                    process.join()
                    process.close()
                    gc.collect()
                    active_processes.remove(process)
                    pbar.update(1)

        active_processes = []
        if self.with_tqdm:
            pbar = tqdm.tqdm(desc=self.pool_name, total=self._max_queue_len)
        else:
            pbar = _DumbTQDM()
        while len(self._job_queue) > 0 and not self._is_terminated:
            while len(self._job_queue) > 0 and len(active_processes) < self.pool_size:
                active_processes.append(self._job_queue.pop(0))
                active_processes[-1].start()
            _scan_through_process()
            time.sleep(0.1)
        while len(active_processes) > 0 and not self._is_terminated:
            _scan_through_process()
            time.sleep(0.1)
        pbar.close()
        self._is_terminated = True

    def __len__(self):
        return len(self._job_queue)

    def stop(self):
        self._is_terminated = True

    def add(self, mp_instance: multiprocessing.Process):
        self._job_queue.append(mp_instance)
        self._max_queue_len += 1

    def append(self, *args, **kwargs):
        self.add(*args, **kwargs)

    def close(self):
        self.join()
