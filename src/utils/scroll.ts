/**
 * 页面滚动控制：回到页首 + 弹层/抽屉的滚动锁定。
 * 两种锁语义不同，不得互相替换：
 * - lockScrollOverflow：隐藏 overflow（移动端 offcanvas 菜单用）
 * - lockScrollPosition：body 钉住 + 记录/恢复滚动位置（右侧抽屉用）
 */

const prefersReducedMotion = (): boolean =>
  typeof window !== 'undefined' &&
  typeof window.matchMedia === 'function' &&
  window.matchMedia('(prefers-reduced-motion: reduce)').matches

/** 回到页首（SSR 环境安全；reduce-motion 时自动改用瞬时滚动） */
export const scrollToTop = (behavior: 'auto' | 'smooth' = 'smooth') => {
  if (typeof window === 'undefined') return
  window.scrollTo({ top: 0, behavior: prefersReducedMotion() ? 'auto' : behavior })
}

/** 当前是否处于滚动锁定状态（body overflow hidden） */
export const isScrollLocked = (): boolean =>
  typeof document !== 'undefined' &&
  !!document.body?.style &&
  document.body.style.overflow === 'hidden'

/**
 * 平滑滚动到目标元素：吸顶 header 偏移感知、reduce-motion 瞬时滚动、
 * 滚动锁定期间延迟执行（OnThisPage/TocDrawer 锚点跳转共用）。
 */
export const scrollToHeading = (
  el: Element,
  offset = 0,
  opts: { lockedDelay?: number } = {},
): void => {
  const top = Math.max(0, el.getBoundingClientRect().top + window.scrollY - offset)
  const doScroll = () => {
    window.scrollTo({ top, behavior: prefersReducedMotion() ? 'auto' : 'smooth' })
  }
  try {
    if (isScrollLocked()) {
      setTimeout(doScroll, opts.lockedDelay ?? 80)
    } else {
      doScroll()
    }
  } catch {
    doScroll()
  }
}

function setDocumentEl(overflow: string, overscrollBehavior: string) {
  const docEl = document.documentElement
  if (docEl) {
    docEl.style.overflow = overflow
    docEl.style.overscrollBehavior = overscrollBehavior
  }
}

function setBody(overflow: string, overscrollBehavior: string) {
  const body = document.body
  if (body) {
    body.style.overflow = overflow
    body.style.overscrollBehavior = overscrollBehavior
  }
}

export function lockScrollOverflow() {
  setDocumentEl('hidden', 'contain')
  setBody('hidden', 'contain')
}

export function unlockScrollOverflow() {
  setDocumentEl('', '')
  setBody('', '')
}

export function lockScrollPosition(): number {
  const scrollY = window.scrollY || window.pageYOffset || document.documentElement.scrollTop || 0
  try {
    const body = document.body
    if (body) {
      body.style.position = 'fixed'
      body.style.top = `-${scrollY}px`
      body.style.left = '0'
      body.style.right = '0'
      body.style.overflow = 'hidden'
    }
    setDocumentEl('', 'contain')
  } catch {
    setDocumentEl('hidden', 'contain')
    setBody('hidden', 'contain')
  }
  return scrollY
}

export function unlockScrollPosition(lockedScrollY: number | null) {
  try {
    const body = document.body
    if (body) {
      body.style.position = ''
      body.style.top = ''
      body.style.left = ''
      body.style.right = ''
      body.style.overflow = ''
    }
    setDocumentEl('', '')
    if (typeof lockedScrollY === 'number') {
      window.scrollTo(0, lockedScrollY)
    }
  } catch {
    setDocumentEl('', '')
    setBody('', '')
  }
}
