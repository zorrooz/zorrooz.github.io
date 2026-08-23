import { ref } from 'vue'
import { defineStore } from 'pinia'
import i18n from '@/i18n/index'

import { HTML_LANG, localeFromPath, type SupportedLocale } from '@/config'
import { currentLocale, persistLocale } from '@/i18n/locale'

export const useLocaleStore = defineStore('locale', () => {
  const locale = ref<SupportedLocale>(currentLocale())

  const setLocale = (nextLocale: SupportedLocale) => {
    locale.value = nextLocale
    i18n.global.locale.value = nextLocale
    persistLocale(nextLocale)
    if (typeof document !== 'undefined') document.documentElement.lang = HTML_LANG[nextLocale]
  }

  const initLocale = () => {
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
