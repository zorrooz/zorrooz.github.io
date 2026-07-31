import fs from 'node:fs'
import path from 'node:path'
import { fileURLToPath, URL } from 'node:url'

import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'
import vueDevTools from 'vite-plugin-vue-devtools'

import { contentDev } from './src/utils/contentDevPlugin.js'

function getArticleRoutes() {
  const categoriesPath = path.resolve('src/content/categories.json')
  if (!fs.existsSync(categoriesPath)) return []
  const data = JSON.parse(fs.readFileSync(categoriesPath, 'utf-8'))
  const routes = []
  if (Array.isArray(data)) {
    for (const section of data) {
      for (const item of section?.items ?? []) {
        for (const art of item?.articles ?? []) routes.push(art.articleUrl)
        for (const cat of item?.categories ?? []) {
          for (const art of cat?.articles ?? []) routes.push(art.articleUrl)
        }
      }
    }
  }
  return routes.filter(Boolean)
}

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
  ssgOptions: {
    includedRoutes: (paths) => [
      ...paths.filter((p) => !p.includes(':')),
      ...getArticleRoutes(),
    ],
    onFinished: () => {
      fs.copyFileSync('dist/index.html', 'dist/404.html')
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
