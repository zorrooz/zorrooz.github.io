export const SITE = {
  url: 'https://zorrooz.github.io',
  author: 'zorrooz',
  email: 'zorrooz@163.com',
  github: 'https://github.com/zorrooz',
  startYear: 2025,
} as const

export const SUPPORTED_LOCALES = ['zh-CN', 'en-US'] as const
export type SupportedLocale = (typeof SUPPORTED_LOCALES)[number]

export const LOCALE_SEGMENTS = ['zh', 'en'] as const
export type LocaleSegment = (typeof LOCALE_SEGMENTS)[number]

export const THEME_MODES = ['auto', 'light', 'dark'] as const

export const DEFAULT_LOCALE: SupportedLocale = 'zh-CN'

/** URL 段（/zh、/en）→ 完整 locale 的唯一映射 */
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

/** Map a URL locale segment (`/zh`, `/en`) or undefined to a supported locale. */
export const localeFromSegment = (segment: string | undefined): SupportedLocale =>
  segment && segment in LOCALE_MAP ? LOCALE_MAP[segment as LocaleSegment] : DEFAULT_LOCALE

/** URL segment for a locale (`/zh` or `/en`). */
export const localeSegmentOf = (locale: SupportedLocale): LocaleSegment => SEGMENT_OF[locale]

/** Extract the locale from a URL path prefix (`/zh/...`, `/en/...`); null when absent. */
export const localeFromPath = (path: string): SupportedLocale | null => {
  const m = path.match(/^\/(zh|en)(?=\/|$)/)
  return m ? LOCALE_MAP[m[1] as LocaleSegment] : null
}

/** Strip a leading `/zh`/`/en` prefix from a path (used by App.vue for SEO links). */
export const stripLocalePrefix = (path: string): string => path.replace(/^\/(zh|en)(?=\/|$)/, '')

/** Coerce any runtime string (i18n locale, storage value) to a supported locale. */
export const toSupportedLocale = (value: string | null | undefined): SupportedLocale =>
  value === SUPPORTED_LOCALES[1] ? SUPPORTED_LOCALES[1] : SUPPORTED_LOCALES[0]

/**
 * Client-preferred locale segment: persisted `locale` takes precedence,
 * then `navigator.language`; falls back to `zh`.
 */
export const preferredLocaleSegment = (): LocaleSegment => {
  if (typeof window !== 'undefined') {
    const saved = localStorage.getItem('locale')
    if (saved === SUPPORTED_LOCALES[1]) return SEGMENT_OF[SUPPORTED_LOCALES[1]]
    if (saved === SUPPORTED_LOCALES[0]) return SEGMENT_OF[SUPPORTED_LOCALES[0]]
    if (navigator.language && navigator.language.toLowerCase().startsWith('en')) return 'en'
  }
  return 'zh'
}
