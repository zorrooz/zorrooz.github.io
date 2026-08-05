/// <reference types="vite/client" />

declare module 'bootstrap'

/** SSR prerender 期间由 main.ts 注入的 locale（contentLoader 在无 router 上下文时读取） */
declare var __GBLOG_LOCALE__: string | undefined
