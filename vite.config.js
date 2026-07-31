import { fileURLToPath, URL } from 'node:url'

import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'
import vueDevTools from 'vite-plugin-vue-devtools'

import { contentDev } from './src/utils/contentDevPlugin.js'

// https://vite.dev/config/
export default defineConfig({
  plugins: [
    vue(),
    contentDev(),
    process.env.NODE_ENV !== 'production' && vueDevTools(),
  ].filter(Boolean),
  server: {
    watch: {
      ignored: ['**/src/content/**'],
    },
  },
  css: {
    preprocessorOptions: {
      scss: {
        quietDeps: true
      }
    }
  },
  resolve: {
    alias: {
      '@': fileURLToPath(new URL('./src', import.meta.url)),
    },
  },
})
