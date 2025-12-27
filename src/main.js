//main.js
import { createApp } from 'vue'
import { createPinia } from 'pinia'
import { RouterLink } from 'vue-router'

import App from './App.vue'
import router from './router'
import 'bootstrap'
import './stores/styles/global.scss'
import i18n from './stores/i18n'
import 'katex/dist/katex.min.css'

import { useAppStore } from './stores/app'

const app = createApp(App)

app.component('RouterLink', RouterLink)

app.use(createPinia())
app.use(router)
app.use(i18n)

const appStore = useAppStore()
appStore.initTheme()
appStore.initLocale()

app.mount('#app')
