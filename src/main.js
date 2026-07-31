//main.js
import { ViteSSG } from 'vite-ssg'
import { createPinia } from 'pinia'
import { VueHeadMixin } from '@unhead/vue'

import App from './App.vue'
import { routes } from './router'
import './stores/styles/global.scss'
import './assets/styles/highlight/github.css'
import './assets/styles/highlight/github-dark-dimmed.css'
import './assets/styles/fa-subset.css'
import i18n from './stores/i18n'
import 'katex/dist/katex.min.css'

import { useAppStore } from './stores/app'

// bootstrap touches document at module scope; SSR build must not include it
if (!import.meta.env.SSR) import('bootstrap')

export const createApp = ViteSSG(
  App,
  { routes },
  ({ app, router, initialState, isClient }) => {
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

    const appStore = useAppStore()
    appStore.initTheme()
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
