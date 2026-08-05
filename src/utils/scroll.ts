/** 回到页首（SSR 环境安全）。Article/PostList/Home 共用的滚动逻辑。 */
export const scrollToTop = (behavior: 'auto' | 'smooth' = 'smooth') => {
  if (typeof window === 'undefined') return
  window.scrollTo({ top: 0, behavior })
}
