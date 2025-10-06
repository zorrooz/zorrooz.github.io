/**
 * 主题管理（最小实现）：
 * - 单一样式入口：src/assets/scss/custom.scss
 * - 切换主题仅设置 document.documentElement 的 data-bs-theme 和持久化到 localStorage
 */

export function setTheme(theme = 'light') {
  const mode = theme === 'dark' ? 'dark' : 'light';
  localStorage.setItem('theme', mode);
  document.documentElement.setAttribute('data-bs-theme', mode);
}

export function initTheme() {
  const saved = localStorage.getItem('theme') || 'light';
  setTheme(saved);
}