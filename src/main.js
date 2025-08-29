import { createApp } from 'vue'
import { createPinia } from 'pinia'
import { RouterLink } from 'vue-router'

import App from './App.vue'
import router from './router'
import 'bootstrap/dist/css/bootstrap.min.css'
import './assets/scss/custom.scss';
import 'bootstrap'
import i18n from './i18n'

const app = createApp(App)

app.component('RouterLink', RouterLink)

app.use(createPinia())
app.use(router)
app.use(i18n)
app.mount('#app')
