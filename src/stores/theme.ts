import { ref } from 'vue'
import { defineStore } from 'pinia'

import { THEME_MODES, THEME_STORAGE_KEY, type ThemeMode } from '@/config'

const setTheme = (mode: ThemeMode) => {
  if (typeof window === 'undefined' || typeof document === 'undefined') return
  const html = document.documentElement

  if (mode === 'auto') {
    const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches
    html.setAttribute('data-bs-theme', prefersDark ? 'dark' : 'light')
    localStorage.removeItem(THEME_STORAGE_KEY)
  } else {
    html.setAttribute('data-bs-theme', mode)
    localStorage.setItem(THEME_STORAGE_KEY, mode)
  }
}

const initialTheme = (): ThemeMode => {
  if (typeof window === 'undefined') return 'auto'
  const saved = localStorage.getItem(THEME_STORAGE_KEY)
  return saved === 'light' || saved === 'dark' ? saved : 'auto'
}

export const useThemeStore = defineStore('theme', () => {
  const theme = ref<ThemeMode>(initialTheme())

  const toggleTheme = () => {
    const currentIndex = THEME_MODES.indexOf(theme.value)
    const nextIndex = (currentIndex + 1) % THEME_MODES.length
    theme.value = THEME_MODES[nextIndex]
    setTheme(theme.value)
  }

  const initTheme = () => {
    setTheme(theme.value)
  }

  return { theme, toggleTheme, initTheme }
})
