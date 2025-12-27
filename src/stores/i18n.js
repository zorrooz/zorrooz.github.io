import { createI18n } from 'vue-i18n'
import zhCN from './locales/zh-CN.js'
import enUS from './locales/en-US.js'

const savedLocale = localStorage.getItem('locale') || 'zh-CN'
document.documentElement.lang = savedLocale

const i18n = createI18n({
  locale: savedLocale,
  fallbackLocale: 'zh-CN',
  messages: {
    'zh-CN': zhCN,
    'en-US': enUS
  },
  legacy: false,
  globalInjection: true
})

export default i18n