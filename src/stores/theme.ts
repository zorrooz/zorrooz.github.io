import { ref } from 'vue'
import { defineStore } from 'pinia'

import { THEME_MODES } from '@/config'

const setTheme = (mode: string) => {
  if (typeof window === 'undefined' || typeof document === 'undefined') return
  const html = document.documentElement

  if (mode === 'auto') {
    const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches
    html.setAttribute('data-bs-theme', prefersDark ? 'dark' : 'light')
    localStorage.removeItem('theme')
  } else {
    html.setAttribute('data-bs-theme', mode)
    localStorage.setItem('theme', mode)
  }
}

export const useThemeStore = defineStore('theme', () => {
  const theme = ref<string>(
    typeof window !== 'undefined' ? localStorage.getItem('theme') || 'auto' : 'auto',
  )

  const toggleTheme = () => {
    const currentIndex = THEME_MODES.indexOf(theme.value as (typeof THEME_MODES)[number])
    const nextIndex = (currentIndex + 1) % THEME_MODES.length
    theme.value = THEME_MODES[nextIndex]
    setTheme(theme.value)
  }

  const initTheme = () => {
    setTheme(theme.value)
  }

  return { theme, toggleTheme, initTheme }
})
