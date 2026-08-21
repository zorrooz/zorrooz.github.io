import { createI18n } from 'vue-i18n'
import zhCN from './locales/zh-CN'
import enUS from './locales/en-US'

import { DEFAULT_LOCALE, HTML_LANG } from '@/config'
import { currentLocale } from '@/locale'

const savedLocale = currentLocale()
if (typeof document !== 'undefined') document.documentElement.lang = HTML_LANG[savedLocale]

const i18n = createI18n({
  locale: savedLocale,
  fallbackLocale: DEFAULT_LOCALE,
  messages: {
    'zh-CN': zhCN,
    'en-US': enUS,
  },
  legacy: false,
  globalInjection: true,
})

export default i18n
