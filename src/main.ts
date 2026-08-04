//main.ts
import { ViteSSG } from 'vite-ssg'
import { createPinia } from 'pinia'
import { VueHeadMixin } from '@unhead/vue'

import App from './App.vue'
import { routes } from './router/index.ts'
import './stores/styles/global.scss'
import './assets/styles/highlight/github.css'
import './assets/styles/highlight/github-dark-dimmed.css'
import '@fortawesome/fontawesome-free/css/fontawesome.min.css'
import '@fortawesome/fontawesome-free/css/solid.min.css'
import '@fortawesome/fontawesome-free/css/brands.min.css'
import i18n from './stores/i18n.ts'
import 'katex/dist/katex.min.css'

import { useAppStore } from './stores/app.ts'
import { reveal } from './utils/reveal'

import { localeFromPrefix } from './config/site.ts'

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

    const appStore = useAppStore()
    appStore.initTheme()
    // SSR: prerender each page in its own locale (/zh or /en path prefix)
    if (import.meta.env.SSR) {
      if (typeof routePath === 'string') {
        const m = routePath.match(/^\/(zh|en)(\/|$)/)
        if (m) {
          const locale = localeFromPrefix(m[1])
          appStore.setLocale(locale)
          // contentLoader (no router context) reads this during prerender
          ;(globalThis as { __GBLOG_LOCALE__?: string }).__GBLOG_LOCALE__ = locale
        }
      }
    }
    // SSR always renders zh-CN; restore the persisted locale on the client
    if (isClient) appStore.initLocale()
    // PWA: register the service worker on every page (index.html carries the manifest link)
    if (isClient && import.meta.env.PROD && 'serviceWorker' in navigator) {
      window.addEventListener('load', () => {
        navigator.serviceWorker.register('/sw.js').catch(() => {})
      })
    }
  },
)
