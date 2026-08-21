export const SITE = {
  url: 'https://zorrooz.github.io',
  author: 'zorrooz',
  email: 'zorrooz@163.com',
  github: 'https://github.com/zorrooz',
  startYear: 2025,
} as const

/** 品牌名（SEO title 后缀，非翻译内容） */
export const BLOG_TITLE = 'zorrooz’s blog'

export const SUPPORTED_LOCALES = ['zh-CN', 'en-US'] as const
export type SupportedLocale = (typeof SUPPORTED_LOCALES)[number]

export const LOCALE_SEGMENTS = ['zh', 'en'] as const
export type LocaleSegment = (typeof LOCALE_SEGMENTS)[number]

export const THEME_MODES = ['auto', 'light', 'dark'] as const
export type ThemeMode = (typeof THEME_MODES)[number]

export const DEFAULT_LOCALE: SupportedLocale = 'zh-CN'

/** localStorage 键：locale（与 html lang / 重定向共用） */
export const LOCALE_STORAGE_KEY = 'locale'
/** localStorage 键：theme */
export const THEME_STORAGE_KEY = 'theme'

/** 英文内容身份后缀（文件与 URL 的一部分；与 scripts/dataConfig.ts 的 EN_SUFFIX 语义一致） */
export const EN_SUFFIX = '-en'

/** 文章路由前缀（路由 / SEO / 生成器共用，改版只改此处） */
export const ARTICLE_ROUTE_PREFIX = '/article'

/** 吸顶 header 高度（锚点滚动偏移，RenderMarkdown/OnThisPage/TocDrawer 共用） */
export const HEADER_OFFSET = 88

/** URL 段（/zh、/en）→ 完整 locale */
export const LOCALE_MAP: Record<LocaleSegment, SupportedLocale> = {
  zh: 'zh-CN',
  en: 'en-US',
}

/** 完整 locale → URL 段 */
export const SEGMENT_OF: Record<SupportedLocale, LocaleSegment> = {
  'zh-CN': 'zh',
  'en-US': 'en',
}

/** html lang 属性值 */
export const HTML_LANG: Record<SupportedLocale, string> = {
  'zh-CN': 'zh-CN',
  'en-US': 'en',
}

/** URL 段（`/zh`/`/en`）或 undefined → 支持的 locale */
export const localeSegmentOf = (locale: SupportedLocale): LocaleSegment => SEGMENT_OF[locale]

/** 从 URL 路径前缀（`/zh/...`、`/en/...`）提取 locale；无前缀返回 null */
export const localeFromPath = (path: string): SupportedLocale | null => {
  const m = path.match(/^\/(zh|en)(?=\/|$)/)
  return m ? LOCALE_MAP[m[1] as LocaleSegment] : null
}

/** 剥离路径前的 `/zh`/`/en` 前缀（App.vue 生成 SEO 链接用） */
export const stripLocalePrefix = (path: string): string => path.replace(/^\/(zh|en)(?=\/|$)/, '')

/** 任意运行时字符串（i18n locale、storage 值）→ 受支持的 locale */
export const toSupportedLocale = (value: string | null | undefined): SupportedLocale =>
  value === SUPPORTED_LOCALES[1] ? SUPPORTED_LOCALES[1] : SUPPORTED_LOCALES[0]

/**
 * 客户端首选 locale 段：优先持久化的 `locale`，其次 `navigator.language`，最后回退 `zh`。
 */
export const preferredLocaleSegment = (): LocaleSegment => {
  if (typeof window !== 'undefined') {
    const saved = localStorage.getItem(LOCALE_STORAGE_KEY)
    if (saved === SUPPORTED_LOCALES[1]) return SEGMENT_OF[SUPPORTED_LOCALES[1]]
    if (saved === SUPPORTED_LOCALES[0]) return SEGMENT_OF[SUPPORTED_LOCALES[0]]
    if (navigator.language && navigator.language.toLowerCase().startsWith('en')) return 'en'
  }
  return 'zh'
}
