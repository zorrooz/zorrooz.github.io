export const SITE = {
  url: 'https://zorrooz.github.io',
  name: 'gblog',
  author: 'zorrooz',
  email: 'zorrooz@163.com',
  github: 'https://github.com/zorrooz',
  startYear: 2025,
} as const

export const SUPPORTED_LOCALES = ['zh-CN', 'en-US'] as const
export type SupportedLocale = (typeof SUPPORTED_LOCALES)[number]

export const LOCALE_PREFIXES = ['/zh', '/en'] as const

export const THEME_MODES = ['auto', 'light', 'dark'] as const

/** Map a URL locale segment (`/zh`, `/en`) or undefined to a supported locale. */
export const localeFromPrefix = (prefix: string | undefined): SupportedLocale => {
  return prefix === 'en' ? SUPPORTED_LOCALES[1] : SUPPORTED_LOCALES[0]
}
