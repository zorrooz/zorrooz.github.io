export const setTheme = (mode) => {
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

export const initTheme = () => {
  const savedMode = localStorage.getItem('theme');
  setTheme(savedMode || 'auto');
};