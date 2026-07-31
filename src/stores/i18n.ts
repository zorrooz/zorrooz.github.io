import { createI18n } from 'vue-i18n'
import zhCN from './locales/zh-CN.ts'
import enUS from './locales/en-US.ts'

import { SUPPORTED_LOCALES } from '@/config/site.ts'

const savedLocale =
  typeof window !== 'undefined'
    ? localStorage.getItem('locale') || SUPPORTED_LOCALES[0]
    : SUPPORTED_LOCALES[0]
if (typeof document !== 'undefined') document.documentElement.lang = savedLocale

const i18n = createI18n({
  locale: savedLocale,
  fallbackLocale: SUPPORTED_LOCALES[0],
  messages: {
    [SUPPORTED_LOCALES[0]]: zhCN,
    [SUPPORTED_LOCALES[1]]: enUS
  },
  legacy: false,
  globalInjection: true
})

export default i18n
