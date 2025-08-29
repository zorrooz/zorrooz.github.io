import { createI18n } from 'vue-i18n'
import zhCN from '../locales/zh-CN.js'
import enUS from '../locales/en-US.js'

// 从 localStorage 读取用户选择，或默认中文
const savedLocale = localStorage.getItem('locale') || 'zh-CN'
document.documentElement.lang = savedLocale // 设置 HTML lang 属性

const i18n = createI18n({
  locale: savedLocale,     // 当前语言
  fallbackLocale: 'zh-CN', // 备用语言
  messages: {
    'zh-CN': zhCN,
    'en-US': enUS
  },
  legacy: false,           // 使用 Composition API 模式
  globalInjection: true    // 允许在模板中使用 $t
})

export default i18n
