/**
 * Sets the theme for the entire application.
 * @param {'light' | 'dark' | 'auto'} mode - The theme mode to set.
 */
export const setTheme = (mode) => {
  const html = document.documentElement;

  if (mode === 'auto') {
    // Follow system preference
    const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
    html.setAttribute('data-bs-theme', prefersDark ? 'dark' : 'light');
    localStorage.removeItem('theme');
  } else {
    // Apply light or dark mode
    html.setAttribute('data-bs-theme', mode);
    localStorage.setItem('theme', mode);
  }
};

/**
 * Initializes the theme based on user preference or system settings.
 */
export const initTheme = () => {
  const savedMode = localStorage.getItem('theme');
  // If there's a saved mode, use it. Otherwise, default to 'auto'.
  setTheme(savedMode || 'auto');
};