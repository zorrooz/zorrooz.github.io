// @ts-check
/** @param {string} mode */
export const setTheme = (mode) => {
  if (typeof window === 'undefined' || typeof document === 'undefined') return
  const html = document.documentElement;

  if (mode === 'auto') {
    const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
    html.setAttribute('data-bs-theme', prefersDark ? 'dark' : 'light');
    localStorage.removeItem('theme');
  } else {
    html.setAttribute('data-bs-theme', mode);
    localStorage.setItem('theme', mode);
  }
};

