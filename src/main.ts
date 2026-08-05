import { createPinia } from 'pinia'
import { VueHeadMixin } from '@unhead/vue'
import { ViteSSG } from 'vite-ssg'

import App from './App.vue'
import { routes } from './router'
import i18n from './i18n'
import { useLocaleStore } from './stores/locale'
import { useThemeStore } from './stores/theme'
import { reveal } from './utils/reveal'
import { localeFromPath } from './config/site'

import './assets/styles/global.scss'
import './assets/styles/highlight/github.css'
import './assets/styles/highlight/github-dark-dimmed.css'
import '@fortawesome/fontawesome-free/css/fontawesome.min.css'
import '@fortawesome/fontawesome-free/css/solid.min.css'
import '@fortawesome/fontawesome-free/css/brands.min.css'
import 'katex/dist/katex.min.css'

// bootstrap touches document at module scope; SSR build must not include it
if (!import.meta.env.SSR) import('bootstrap')

export const createApp = ViteSSG(
  App,
  { routes },
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
    // SSR: prerender each page in its own locale (/zh or /en path prefix)
    if (import.meta.env.SSR && typeof routePath === 'string') {
      const locale = localeFromPath(routePath)
      if (locale) {
        localeStore.setLocale(locale)
        // contentLoader (no router context) reads this during prerender
        globalThis.__GBLOG_LOCALE__ = locale
      }
    }
    if (isClient) localeStore.initLocale()
    // PWA: register the service worker on every page (index.html carries the manifest link)
    if (isClient && import.meta.env.PROD && 'serviceWorker' in navigator) {
      window.addEventListener('load', () => {
        navigator.serviceWorker.register('/sw.js').catch(() => {})
      })
    }
  },
)
