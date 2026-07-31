import i18n from '@/stores/i18n.ts'

import { LOCALE_PREFIXES, SUPPORTED_LOCALES } from '@/config/site.ts'

/**
 * Prefix an internal path with the current locale (/zh or /en).
 * Accepts '/category' style paths or bare article paths like
 * 'notes/Omics/genomics/bwa/bwa' (the leading slash is added).
 */
export const toLocalePath = (path: string): string => {
  const prefix = i18n.global.locale.value === SUPPORTED_LOCALES[1] ? LOCALE_PREFIXES[1] : LOCALE_PREFIXES[0]
  return `${prefix}${path.startsWith('/') ? path : `/${path}`}`
}
