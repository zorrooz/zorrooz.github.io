import fs from 'node:fs'
import path from 'node:path'
import { fileURLToPath, URL } from 'node:url'

import type { UserConfig } from 'vite'
import vue from '@vitejs/plugin-vue'
import vueDevTools from 'vite-plugin-vue-devtools'
import { VitePWA } from 'vite-plugin-pwa'
import type { ViteSSGOptions } from 'vite-ssg'

import { contentDev } from './vite/contentDevPlugin.ts'
import { dataDir, contentDir } from './scripts/dataConfig.ts'

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
    const categoriesPath = path.resolve(contentDir, file)
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
    ...collect('categories.json').map((r) => `/zh${r}`),
    ...collect('categories-en.json').map((r) => `/en${r}`),
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
      // 数据目录在仓库外（git worktree ../blog-data），由 contentDevPlugin 负责监听重生成，
      // 这里忽略避免与 Vite 自身 watcher 双重重载冲突
      ignored: [`${dataDir}/**`],
    },
    fs: {
      // 允许 dev server 读取仓库外数据目录（@data 别名指向它）
      allow: [dataDir, fileURLToPath(new URL('.', import.meta.url))],
    },
  },
  ssgOptions: {
    // i18n.global.locale 是模块单例，串行渲染保证各页 locale 不被并发页覆盖
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
      // 数据分支目录：所有 content*/contentHtml 通过 @data 引用（见 contentLoader）
      '@data': dataDir,
      // vue-i18n v12 alpha 的 browser 条件会解析到 runtime-only 版（无 message compiler），
      // 导致生产构建中 t() 不编译插值占位符（页面显示 {count}/{minutes}）。
      // 强制使用完整版（带运行时编译器），与 dev/SSG 行为一致。
      'vue-i18n': fileURLToPath(
        new URL('./node_modules/vue-i18n/dist/vue-i18n.esm-browser.js', import.meta.url),
      ),
    },
  },
}

export default config
