/**
 * 「当前 locale」的单一事实来源。
 * 取值顺序：SSR/路由注入的 `__GBLOG_LOCALE__`（main.ts 设置）→ URL 前缀 → localStorage → DEFAULT_LOCALE。
 * URL 优先于 localStorage：直接打开 /zh/ 链接时内容必须与 URL 一致，
 * 否则 UI 按 URL、数据按存储各走一套，产生 /zh/...-en 之类的混合链接。
 * 所有消费方（contentLoader / i18n 初始化 / locale store）一律走这里，避免各自实现漂移。
 */
import {
  DEFAULT_LOCALE,
  LOCALE_STORAGE_KEY,
  SUPPORTED_LOCALES,
  localeFromPath,
  toSupportedLocale,
  type SupportedLocale,
} from '@/config'

const injectedLocale = (): SupportedLocale | null => {
  const v = (globalThis as { __GBLOG_LOCALE__?: unknown }).__GBLOG_LOCALE__
  return v === SUPPORTED_LOCALES[0] || v === SUPPORTED_LOCALES[1] ? v : null
}

/** 当前生效 locale（SSR 安全） */
export const currentLocale = (): SupportedLocale => {
  if (typeof window !== 'undefined') {
    return (
      injectedLocale() ??
      localeFromPath(window.location.pathname) ??
      toSupportedLocale(localStorage.getItem(LOCALE_STORAGE_KEY))
    )
  }
  return DEFAULT_LOCALE
}

/** 持久化 locale（store.setLocale 调用；写 localStorage 并同步 html lang 由调用方负责） */
export const persistLocale = (locale: SupportedLocale): void => {
  if (typeof window !== 'undefined') localStorage.setItem(LOCALE_STORAGE_KEY, locale)
}
