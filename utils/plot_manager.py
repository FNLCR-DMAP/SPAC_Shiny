from shiny import ui, reactive
import queue as queue_module
import threading
import multiprocessing


class PlotManager:
    """
    Manages a single plot's lifecycle: starting, monitoring, and cancelling.
    Reduces boilerplate in each server file to just the plot-specific logic.
    Writes to shared['plot_registry'] so app.py cancel_plot() can terminate it.
    """

    def __init__(self, plot_id, shared, plot_type='process', data_key=None):
        """
        Parameters
        ----------
        plot_id : str
            Unique ID matching the STOP_BUTTON_MAP key in app.py
            e.g. 'boxplot', 'spatial', 'histogram'
        shared : dict
            The shared reactive values dictionary from app.py
        plot_type : str
            'process' for multiprocessing (matplotlib plots)
            'thread' for threading (widget/Plotly plots)
        data_key : str or None
            Key in shared to clear on cancel e.g. 'df_boxplot', or None
        """
        self.plot_id = plot_id
        self.shared = shared
        self.plot_type = plot_type
        self.data_key = data_key

        self.current_proc = reactive.Value(None)
        self.result = reactive.Value(None)
        self.is_calculating = reactive.Value(False)
        self.result_queue = None

    def _register(self, proc, q):
        """Register this plot in the global registry for app.py to cancel."""
        self.shared['plot_registry'][self.plot_id] = {
            'type': self.plot_type,
            'proc': proc,
            'queue': q,
            'is_calculating': self.is_calculating,
            'result': self.result,
            'current_proc_reactive': self.current_proc,
            'data_key': self.data_key
        }

    def _unregister(self):
        """Remove this plot from the global registry when done."""
        self.shared['plot_registry'].pop(self.plot_id, None)

    def start_thread(self, target_fn):
        """
        Start a threading.Thread for widget/Plotly plots.
        target_fn should use self.result_queue to return results.

        Parameters
        ----------
        target_fn : callable
            Function with no args that puts result into self.result_queue
        """
        self.result.set(None)
        if self.data_key:
            self.shared[self.data_key].set(None)

        self.result_queue = queue_module.Queue()
        self.current_proc.set(None)

        thread = threading.Thread(target=target_fn, daemon=True)
        thread.start()
        self.current_proc.set(thread)
        self._register(thread, self.result_queue)
        self.is_calculating.set(True)

    def start_process(self, target_fn, args=()):
        """
        Start a multiprocessing.Process for matplotlib plots.
        target_fn signature must be: target_fn(queue, *args)

        Parameters
        ----------
        target_fn : callable
            Worker function - queue is automatically prepended to args
        args : tuple
            Additional arguments passed after queue
        """
        self.result.set(None)
        if self.data_key:
            self.shared[self.data_key].set(None)

        q = multiprocessing.Queue()
        self.result_queue = q
        self.current_proc.set(None)

        p = multiprocessing.Process(target=target_fn, args=(q,) + args)
        p.start()
        self.current_proc.set((p, q))
        self._register(p, q)
        self.is_calculating.set(True)

    def check_thread(self, on_result):
        """
        Poll thread status. Call this inside a reactive Effect with no event.
        Uses reactive.invalidate_later(1.0) to keep polling.

        Parameters
        ----------
        on_result : callable
            Called with the result when thread completes successfully
        """
        thread_info = self.current_proc.get()
        if thread_info is None:
            return

        if not thread_info.is_alive() or not self.result_queue.empty():
            try:
                res = self.result_queue.get_nowait()
                if isinstance(res, str) and res.startswith("Error"):
                    ui.notification_show(res, type="error")
                else:
                    on_result(res)
            except queue_module.Empty:
                pass
            self.current_proc.set(None)
            self._unregister()
            self.is_calculating.set(False)
        else:
            reactive.invalidate_later(1.0)

    def check_process(self, on_result):
        """
        Poll process status. Call this inside a reactive Effect with no event.
        Uses reactive.invalidate_later(1.0) to keep polling.

        Parameters
        ----------
        on_result : callable
            Called with the result when process completes successfully
        """
        p_info = self.current_proc.get()
        if p_info is None:
            return

        p, q = p_info
        if not p.is_alive() or not q.empty():
            try:
                res = q.get(timeout=1)
                if isinstance(res, str) and res.startswith("Error"):
                    ui.notification_show(res, type="error")
                else:
                    on_result(res)
            except Exception:
                pass
            p.join()
            self.current_proc.set(None)
            self._unregister()
            self.is_calculating.set(False)
        else:
            reactive.invalidate_later(1.0)

    def stop_button_ui(self, button_id, label="Cancel Render",
                       class_="btn-danger", style="width: 180px;", **kwargs):
        """
        Return a cancel button if calculating, else None.
        Use inside a @render.ui function.

        Parameters
        ----------
        button_id : str
            The input ID for the stop button e.g. 'stop_boxplot'
        label : str
            Button label text
        class_ : str
            CSS class for button styling
        style : str
            Inline CSS styling
        **kwargs
            Additional HTML attributes to pass to button
        """
        if self.is_calculating.get():
            return ui.input_action_button(
                button_id, label,
                class_=class_, style=style, **kwargs
            )
        return None

    def download_button_ui(self, download_id, label="Download Data",
                           class_="btn-warning", style="width: 180px;", **kwargs):
        """
        Return a download button if data exists and not calculating, else None.
        Use inside a @render.ui function.

        Parameters
        ----------
        download_id : str
            The input ID for the download button e.g. 'download_boxplot'
        label : str
            Button label text
        class_ : str
            CSS class for button styling
        style : str
            Inline CSS styling
        **kwargs
            Additional HTML attributes to pass to button
        """
        if not self.is_calculating.get():
            try:
                if (self.data_key and
                        self.shared.get(self.data_key) and
                        self.shared[self.data_key].get() is not None):
                    return ui.download_button(
                        download_id, label,
                        class_=class_, style=style, **kwargs
                    )
            except (KeyError, AttributeError):
                pass
        return None

    def plot_download_button_ui(self, download_id, label="Download Plot",
                                class_="btn-primary", style="width: 180px;", **kwargs):
        """
        Show download plot button when not calculating and PNG result exists.
        Only works for process-based plots (matplotlib -> PNG bytes).
        Use inside a @render.ui function.

        Parameters
        ----------
        download_id : str
            The input ID for the download button e.g. 'download_histogram2_plot'
        label : str
            Button label text
        class_ : str
            CSS class for button styling
        style : str
            Inline CSS styling
        **kwargs
            Additional HTML attributes to pass to button
        """
        if not self.is_calculating.get():
            try:
                result = self.result.get()
                # Check if result is PNG bytes (not a tuple like boxplot)
                if result is not None and isinstance(result, bytes):
                    return ui.download_button(
                        download_id, label,
                        class_=class_, style=style, **kwargs
                    )
            except (KeyError, AttributeError):
                pass
        return None

    def create_plot_download_handler(self):
        """
        Returns download handler for PNG plot (process-based only).
        The PNG bytes are already in pm.result from the worker.
        Use with @render.download decorator.

        Returns
        -------
        callable
            Function that returns (bytes, mimetype) for download

        Example
        -------
        @render.download(filename="my_plot.png")
        def download_my_plot():
            return pm.create_plot_download_handler()()
        """
        def handler():
            result = self.result.get()
            if result is not None and isinstance(result, bytes):
                return result, "image/png"
            return None
        return handler