import { ref } from 'vue'
import { defineStore } from 'pinia'
import { setTheme } from './theme.ts'
import i18n from './i18n.ts'

import { localeFromPrefix, SUPPORTED_LOCALES, THEME_MODES, type SupportedLocale } from '@/config/site.ts'

export const useAppStore = defineStore('app', () => {
  // 主题状态
  const theme = ref<string>(
    typeof window !== 'undefined' ? localStorage.getItem('theme') || 'auto' : 'auto'
  )

  // 语言状态
  const locale = ref<SupportedLocale>(
    typeof window !== 'undefined'
      ? (localStorage.getItem('locale') as SupportedLocale) || SUPPORTED_LOCALES[0]
      : SUPPORTED_LOCALES[0]
  )

  // 主题切换
  const toggleTheme = () => {
    const currentIndex = THEME_MODES.indexOf(theme.value as (typeof THEME_MODES)[number])
    const nextIndex = (currentIndex + 1) % THEME_MODES.length
    theme.value = THEME_MODES[nextIndex]
    setTheme(theme.value)
  }

  // 语言切换（仅状态；URL 前缀的切换由 AppHeader 导航完成）
  const setLocale = (nextLocale: SupportedLocale) => {
    locale.value = nextLocale

    // 更新i18n和localStorage
    i18n.global.locale.value = nextLocale
    if (typeof window !== 'undefined') localStorage.setItem('locale', nextLocale)
    if (typeof document !== 'undefined') document.documentElement.lang = nextLocale
  }

  // 初始化主题
  const initTheme = () => {
    setTheme(theme.value)
  }

  // 初始化语言：优先取 URL 前缀（/zh /en），回退到 localStorage
  const initLocale = () => {
    if (typeof window !== 'undefined') {
      const m = window.location.pathname.match(/^\/(zh|en)(\/|$)/)
      if (m) {
        locale.value = localeFromPrefix(m[1])
      }
    }
    i18n.global.locale.value = locale.value
    if (typeof document !== 'undefined') document.documentElement.lang = locale.value
  }

  return {
    theme,
    locale,
    setLocale,
    toggleTheme,
    initTheme,
    initLocale
  }
})
