// @ts-check
import i18n from '@/stores/i18n'

/**
 * Prefix an internal path with the current locale (/zh or /en).
 * Accepts '/category' style paths or bare article paths like
 * 'notes/Omics/genomics/bwa/bwa' (the leading slash is added).
 * @param {string} path
 * @returns {string}
 */
export const toLocalePath = (path) => {
  const prefix = i18n.global.locale.value === 'en-US' ? '/en' : '/zh'
  return `${prefix}${path.startsWith('/') ? path : `/${path}`}`
}
