import { createPinia } from 'pinia'
import { VueHeadMixin } from '@unhead/vue'
import { ViteSSG } from 'vite-ssg'

import App from './App.vue'
import { routes, scrollBehavior } from './router'
import i18n from './i18n'
import { reveal } from './directives/reveal'
import { useLocaleStore } from './stores/locale'
import { useThemeStore } from './stores/theme'
import { localeFromPath } from './config'

import './assets/styles/global.scss'
import './assets/styles/highlight/github.css'
import './assets/styles/highlight/github-dark-dimmed.css'
import '@fortawesome/fontawesome-free/css/fontawesome.min.css'
import '@fortawesome/fontawesome-free/css/solid.min.css'
import '@fortawesome/fontawesome-free/css/brands.min.css'
import 'katex/dist/katex.min.css'

/** bootstrap 在模块顶层触碰 document，SSR 构建不得引入 */
if (!import.meta.env.SSR) import('bootstrap')

export const createApp = ViteSSG(
  App,
  { routes, scrollBehavior },
  ({ app, router, initialState, isClient, routePath }) => {
    const pinia = createPinia()
    app.use(pinia)
    if (import.meta.env.SSR) {
      initialState.pinia = pinia.state.value
    } else if (initialState.pinia) {
      pinia.state.value = initialState.pinia
    }

    app.use(router)
    app.use(i18n)
    app.mixin(VueHeadMixin)
    app.directive('reveal', reveal)

    const themeStore = useThemeStore()
    const localeStore = useLocaleStore()
    themeStore.initTheme()
    /** 按 /zh 或 /en 路径前缀逐页预渲染，注入当前 locale（contentLoader 预渲染期读取） */
    if (import.meta.env.SSR && typeof routePath === 'string') {
      const locale = localeFromPath(routePath)
      if (locale) {
        localeStore.setLocale(locale)
        globalThis.__GBLOG_LOCALE__ = locale
      }
    }
    if (isClient) localeStore.initLocale()
    if (isClient && import.meta.env.PROD && 'serviceWorker' in navigator) {
      window.addEventListener('load', () => {
        navigator.serviceWorker.register('/sw.js').catch(() => {})
      })
    }
  },
)
