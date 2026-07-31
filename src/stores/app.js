// @ts-check
import { ref } from 'vue'
import { defineStore } from 'pinia'
import { setTheme } from './theme'
import i18n from './i18n'

export const useAppStore = defineStore('app', () => {
  // 主题状态
  /** @type {import('vue').Ref<string>} */
  const theme = ref(
    typeof window !== 'undefined' ? localStorage.getItem('theme') || 'auto' : 'auto'
  )
  
  // 语言状态
  /** @type {import('vue').Ref<string>} */
  const locale = ref(
    typeof window !== 'undefined' ? localStorage.getItem('locale') || 'zh-CN' : 'zh-CN'
  )
  
  // 主题切换
  const toggleTheme = () => {
    const modes = ['auto', 'light', 'dark']
    const currentIndex = modes.indexOf(theme.value)
    const nextIndex = (currentIndex + 1) % modes.length
    theme.value = modes[nextIndex]
    setTheme(theme.value)
  }
  
  // 语言切换
  const toggleLanguage = () => {
    const locales = ['zh-CN', 'en-US']
    const currentIndex = locales.indexOf(locale.value)
    const nextIndex = (currentIndex + 1) % locales.length
    locale.value = locales[nextIndex]
    
    // 更新i18n和localStorage
    i18n.global.locale.value = /** @type {'zh-CN' | 'en-US'} */ (locale.value)
    localStorage.setItem('locale', locale.value)
    document.documentElement.lang = locale.value
  }
  
  // 初始化主题
  const initTheme = () => {
    setTheme(theme.value)
  }
  
  // 初始化语言
  const initLocale = () => {
    i18n.global.locale.value = /** @type {'zh-CN' | 'en-US'} */ (locale.value)
    if (typeof document !== 'undefined') document.documentElement.lang = locale.value
  }
  
  return {
    theme,
    locale,
    toggleTheme,
    toggleLanguage,
    initTheme,
    initLocale
  }
})