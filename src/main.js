//main.js
import { createApp } from 'vue'
import { createPinia } from 'pinia'
import { RouterLink } from 'vue-router'

import App from './App.vue'
import router from './router'
import 'bootstrap'
import './assets/styles/global.scss'
import i18n from './i18n'
import 'katex/dist/katex.min.css'

import { initTheme } from './utils/theme'

const app = createApp(App)

app.component('RouterLink', RouterLink)

app.use(createPinia())
app.use(router)
app.use(i18n)

/* 初始化主题（根据 localStorage 加载对应 CSS 文件） */
initTheme()

app.mount('#app')
