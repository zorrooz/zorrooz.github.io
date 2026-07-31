import fs from 'node:fs'
import path from 'node:path'
import { fileURLToPath, URL } from 'node:url'

import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'
import vueDevTools from 'vite-plugin-vue-devtools'
import { VitePWA } from 'vite-plugin-pwa'

import { contentDev } from './src/utils/contentDevPlugin.js'

function getArticleRoutes() {
  const collect = (file) => {
    const categoriesPath = path.resolve(file)
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
  return [
    ...collect('src/content/categories.json').map((r) => `/zh${r}`),
    ...collect('src/content/categories-en.json').map((r) => `/en${r}`),
  ]
}

// https://vite.dev/config/
export default defineConfig({
  plugins: [
    vue(),
    contentDev(),
    VitePWA({
      registerType: 'autoUpdate',
      injectRegister: false,
      includeAssets: ['favicon.ico', 'icons/gblog-512.svg'],
      manifest: {
        name: 'gblog',
        short_name: 'gblog',
        description: '个人技术博客：生物信息学、编程与科研笔记',
        lang: 'zh-CN',
        display: 'standalone',
        start_url: '/',
        theme_color: '#047aff',
        background_color: '#ffffff',
        icons: [
          { src: 'favicon.ico', sizes: '48x48', type: 'image/x-icon' },
          { src: 'icons/gblog-512.svg', sizes: '512x512', type: 'image/svg+xml' },
          { src: 'icons/gblog-512.svg', sizes: '512x512', type: 'image/svg+xml', purpose: 'maskable' },
        ],
      },
      workbox: {
        globPatterns: ['**/*.{js,css,woff2,ttf,otf,png,svg,ico,json}'],
        navigateFallback: '/index.html',
        navigateFallbackDenylist: [/^\/sw\.js$/],
        maximumFileSizeToCacheInBytes: 20 * 1024 * 1024,
      },
    }),
    process.env.NODE_ENV !== 'production' && vueDevTools(),
  ].filter(Boolean),
  server: {
    watch: {
      ignored: ['**/src/content/**'],
    },
  },
  ssgOptions: {
    // i18n.global.locale is a module singleton; serial rendering keeps each
    // page's locale from being clobbered by a concurrent page's setup
    concurrency: 1,
    includedRoutes: (paths) => {
      const redirectOnly = new Set(['/', '/category', '/resource', '/about'])
      const prefixed = paths.filter((p) => !p.includes(':') && !redirectOnly.has(p))
      return [
        ...prefixed,
        ...getArticleRoutes(),
      ]
    },
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
