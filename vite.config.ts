import fs from 'node:fs'
import path from 'node:path'
import { fileURLToPath, URL } from 'node:url'

import type { UserConfig } from 'vite'
import vue from '@vitejs/plugin-vue'
import vueDevTools from 'vite-plugin-vue-devtools'
import { VitePWA } from 'vite-plugin-pwa'
import type { ViteSSGOptions } from 'vite-ssg'

import { contentDev } from './src/utils/contentDevPlugin.ts'

interface ArticleEntry {
  articleUrl?: unknown
}

interface CategoryItem {
  articles?: ArticleEntry[]
  categories?: Array<{ articles?: ArticleEntry[] }>
}

interface CategorySection {
  items?: CategoryItem[]
}

function getArticleRoutes(): string[] {
  const collect = (file: string): string[] => {
    const categoriesPath = path.resolve(import.meta.dirname, file)
    if (!fs.existsSync(categoriesPath)) return []
    const data = JSON.parse(fs.readFileSync(categoriesPath, 'utf-8')) as CategorySection[]
    const routes: string[] = []
    if (Array.isArray(data)) {
      for (const section of data) {
        for (const item of section?.items ?? []) {
          for (const art of item?.articles ?? []) routes.push(String(art.articleUrl))
          for (const cat of item?.categories ?? []) {
            for (const art of cat?.articles ?? []) routes.push(String(art.articleUrl))
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
const config: UserConfig & { ssgOptions?: ViteSSGOptions } = {
  plugins: [
    vue(),
    contentDev(),
    VitePWA({
      registerType: 'autoUpdate',
      injectRegister: false,
      includeAssets: ['favicon.png', 'apple-touch-icon.png', 'icon-512.png'],
      manifest: {
        name: 'zorrooz’s blog',
        short_name: 'zorrooz',
        description: '个人技术博客：生物信息学、编程与结构生物学，记录编程教程与科研笔记',
        lang: 'zh-CN',
        display: 'standalone',
        start_url: '/',
        theme_color: '#5aa0f8',
        background_color: '#ffffff',
        icons: [
          { src: '/favicon.png', sizes: '64x64', type: 'image/png' },
          { src: '/apple-touch-icon.png', sizes: '180x180', type: 'image/png' },
          { src: '/icon-512.png', sizes: '512x512', type: 'image/png' },
          { src: '/icon-512.png', sizes: '512x512', type: 'image/png', purpose: 'maskable' },
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
    includedRoutes: (paths: string[]) => {
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
}

export default config
