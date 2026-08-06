import { ref } from 'vue'
import { defineStore } from 'pinia'
import i18n from '@/i18n/index'

import {
  DEFAULT_LOCALE,
  HTML_LANG,
  localeFromPath,
  toSupportedLocale,
  type SupportedLocale,
} from '@/config'

export const useLocaleStore = defineStore('locale', () => {
  const locale = ref<SupportedLocale>(
    typeof window !== 'undefined'
      ? toSupportedLocale(localStorage.getItem('locale'))
      : DEFAULT_LOCALE,
  )

  const setLocale = (nextLocale: SupportedLocale) => {
    locale.value = nextLocale
    i18n.global.locale.value = nextLocale
    if (typeof window !== 'undefined') localStorage.setItem('locale', nextLocale)
    if (typeof document !== 'undefined') document.documentElement.lang = HTML_LANG[nextLocale]
  }

  const initLocale = () => {
    // 优先取 URL 前缀（/zh /en），回退到 localStorage
    if (typeof window !== 'undefined') {
      const localeFromUrl = localeFromPath(window.location.pathname)
      if (localeFromUrl) {
        locale.value = localeFromUrl
      }
    }
    i18n.global.locale.value = locale.value
    if (typeof document !== 'undefined') document.documentElement.lang = HTML_LANG[locale.value]
  }

  return {
    locale,
    setLocale,
    initLocale,
  }
})
