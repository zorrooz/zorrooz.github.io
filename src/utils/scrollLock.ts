/**
 * 弹层/抽屉的页面滚动锁定（两种既有实现，语义不同，不得互相替换）：
 * - lockScrollOverflow：隐藏 overflow（滚动位置天然保持；移动端 offcanvas 菜单用）
 * - lockScrollPosition：body 钉住 + 记录/恢复滚动位置（TocDrawer 右侧抽屉用）
 */

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
